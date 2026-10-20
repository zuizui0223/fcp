"""Biologically consequential transformations in the molecular metadata bridge."""
import importlib.util
import gzip
import math
from pathlib import Path
import tempfile
import unittest
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('expansion', ROOT/'scripts/analysis/build_fcp_molecular_expansion_20261019.py')
MOD=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)

class MolecularExpansionTests(unittest.TestCase):
    def test_names_do_not_collapse_subspecies_or_hybrids(self):
        self.assertEqual(MOD.key('  Matthiola   INCANA '),'matthiola incana')
        self.assertNotEqual(MOD.key('Raphanus sativus var. sativus'),MOD.key('Raphanus sativus'))
        self.assertNotEqual(MOD.key('Rosa hybrid cultivar'),MOD.key('Rosa rugosa'))
        self.assertNotEqual(MOD.key('Verbena'),MOD.key('Verbena stricta'))

    def test_geo_uses_phenotype_attribute_not_white_title_initial(self):
        text='''^SAMPLE = GSM1
!Sample_title = Orychophragmus_violaceus_W1
!Sample_organism_ch1 = Orychophragmus violaceus
!Sample_characteristics_ch1 = phenotype: white-violet flower
!Sample_characteristics_ch1 = developmental stage: half-open flowers
!Sample_relation = BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN00001
!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX00001
'''
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'sample.gz'
            path.write_bytes(gzip.compress(text.encode()))
            frame=MOD.geo_samples(path)
        self.assertEqual(frame.phenotype.iloc[0],'white-violet flower')
        self.assertEqual(frame.sample_accession.iloc[0],'SAMN00001')
        self.assertEqual(frame.experiment_accession.iloc[0],'SRX00001')

    def test_frozen_geo_retains_nonobvious_colors(self):
        frame=MOD.geo_samples(MOD.EVID/'geo_GSE93879_family.soft.gz')
        self.assertEqual(len(frame),44)
        self.assertEqual(frame.loc[frame.sample_title.eq('Matthiola_incana_RB1'),'phenotype'].iloc[0],'pink flower')
        self.assertEqual(frame.loc[frame.sample_title.eq('Matthiola_incana_LP1'),'phenotype'].iloc[0],'modena flower')

    def test_query_union_deduplicates_runs_not_species(self):
        old=MOD.read_runs(MOD.OLD/'ena_title_screen_runs.tsv')
        new=MOD.read_runs(MOD.EVID/'ena_expanded_title_screen_runs.tsv')
        old_ids,new_ids=set(old.run_accession),set(new.run_accession)
        self.assertEqual(len(old_ids & new_ids),12)
        self.assertEqual(len(old_ids | new_ids),801)
        self.assertEqual(new.study_accession.nunique(),65)

    def test_fpkm_identity_limits_and_descriptive_arithmetic(self):
        geo=MOD.geo_samples(MOD.EVID/'geo_GSE93879_family.soft.gz')
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)
            result=MOD.fpkm_descriptives(MOD.EVID,geo,output)
            mapping=pd.read_csv(output/'processed_expression_column_map.csv').fillna('')
            ratios=pd.read_csv(output/'pathway_fpkm_label_conditional_ratios.csv')
            check=pd.read_csv(output/'processed_expression_validation.csv')
        white=mapping.loc[mapping.literal_source_header.eq('Matthiola_incana_W1')]
        self.assertEqual(white.source_header_occurrence.tolist(),[1,2])
        self.assertEqual(white.sample_accession.tolist(),['',''])
        self.assertFalse((mapping.literal_source_header=='Matthiola_incana_W2').any())
        self.assertEqual(check.loc[check.species.eq('Lobularia maritima'),'status'].iloc[0],'blocked_missing_header_separator')
        ans=ratios.loc[ratios.species.eq('Orychophragmus violaceus') & ratios.gene_label.eq('ANS_LDOX_TT18')].iloc[0]
        self.assertEqual(ans.white_label,'white-violet flower')
        self.assertAlmostEqual(ans.log2_white_over_colored_numeric_means, math.log2(((25.6058+67.6249)/2)/((458.167+93.6805)/2)),places=10)
        self.assertEqual(result['count_DE_tests'],0)

if __name__=='__main__':
    unittest.main()
