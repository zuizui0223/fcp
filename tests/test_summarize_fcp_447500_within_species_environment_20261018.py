"""Spatial sampling, complete-case eligibility and within-species effect definitions."""
import unittest
import tempfile
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis import summarize_fcp_447500_within_species_environment_20261018 as m


class EnvironmentDescriptionTests(unittest.TestCase):
    def test_extracted_raster_must_match_authenticated_archive_member(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); raster=root/'climate.tif'; archive=root/'frozen.zip'
            raster.write_bytes(b'original climate raster')
            with zipfile.ZipFile(archive,'w') as z: z.writestr(raster.name,raster.read_bytes())
            m.verify_raster_member(archive,raster)
            raster.write_bytes(b'modified climate raster')
            with self.assertRaisesRegex(ValueError,'archive member'):
                m.verify_raster_member(archive,raster)

    def test_floor_pixels_and_missing_bounds_never_clip(self):
        a = np.array([[10., 20., 30.], [40., -999., 60.]])
        lat = np.array([89.9, 89., 89.1, 90., 88., np.nan, 90.1])
        lon = np.array([-179.9, -179., -177.1, -180., -179., -179., -179.])
        got = m.sample_array(a, lat, lon, -999., step=1.)
        np.testing.assert_allclose(got[:4], [10., np.nan, 30., 10.], equal_nan=True)
        self.assertTrue(np.isnan(got[4:]).all())

    def test_geotiff_requires_epsg4326_pixel_area_and_frozen_transform(self):
        tags = {33550:(1/6,1/6,0),33922:(0,0,0,-180,90,0),
                34735:(1,1,0,3,1024,0,1,2,1025,0,1,1,2048,0,1,4326),42113:'-3.4e38'}
        m.validate_tags(tags, (2160,1080))
        for key, value in [(33550,(1,1,0)),(33922,(0,0,0,0,0,0)),
                           (34735,(1,1,0,1,2048,0,1,3857))]:
            bad = dict(tags); bad[key] = value
            with self.assertRaises(ValueError): m.validate_tags(bad,(2160,1080))

    def test_difference_is_white_minus_other_over_all_colour_population_sd(self):
        x = np.array([1.,3.,5.,7.]); white = np.array([True,True,False,False])
        raw, sd, effect, state = m.effect(x,white)
        self.assertEqual(raw,-4.)
        self.assertAlmostEqual(sd,np.sqrt(5.))
        self.assertAlmostEqual(effect,-4/np.sqrt(5.))
        self.assertEqual(state,'estimable')
        raw, sd, effect, state = m.effect(np.ones(4),white)
        self.assertEqual(raw,0.)
        self.assertTrue(np.isnan(effect))
        self.assertEqual(state,'zero')

    def test_complete_case_gate_is_shared_across_all_variables(self):
        photo = pd.DataFrame({'inat_taxon_id':[1]*42,'species':['A b']*42,
                              'fcp_source_cohort':['discovery']*42,'morph':['white']*6+['red_pink']*36})
        env = np.tile(np.arange(42.)[:,None],(1,len(m.VARIABLES)))
        env[:2,1] = np.nan  # Loss in one variable removes these photos for every variable.
        records, panel = m.species_results(photo,env,minimum_each=5)
        self.assertEqual(len(panel),0)  # Complete cases have only four white photographs.
        records, panel = m.species_results(photo,env,minimum_each=4)
        self.assertEqual(len(panel),1)
        self.assertEqual(len(records),len(m.VARIABLES))
        self.assertTrue((records.n_complete_classifiable==40).all())
        self.assertTrue((records.n_white==4).all())

    def test_summary_counts_only_estimable_effects_and_is_species_equal(self):
        rows=pd.DataFrame({'source_cohort':['discovery']*3,'variable':['BIO1']*3,
                           'standardized_difference':[1.,-1.,np.nan],
                           'raw_difference':[2.,-2.,0.]})
        summary=m.summarize(rows,'primary')
        q=summary.loc[(summary.group=='pooled') & (summary.metric=='standardized_difference')].iloc[0]
        self.assertEqual(q.n_panel_species,3)
        self.assertEqual(q.n_estimable,2)
        self.assertEqual(q.positive_fraction,.5)
        self.assertEqual(q.zero_fraction,0.)
        self.assertEqual(q['mean'],0.)


if __name__ == '__main__': unittest.main()
