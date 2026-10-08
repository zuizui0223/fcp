"""End-to-end source fixtures and fail-closed tampering for FCP quality queue."""
from __future__ import annotations
import copy
import os
import unittest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
from build_fcp_validation_photo_quality_queue_20261009 import assemble, checked


class QualityQueueTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(os.environ["FCP_QUALITY_SOURCE_ROOT"])
        cls.gap1 = checked(root/"gap1.json", "gap1")
        cls.gap2 = checked(root/"gap2.json", "gap2")
        cls.status = checked(root/"gap2_status.json", "status")
        cls.queue = checked(root/"queue.csv", "queue")

    def test_all_metadata_preserved_and_unreviewed(self):
        r = assemble(self.gap1, self.gap2, self.status, self.queue)
        self.assertEqual((r["n_candidate_species"],r["n_photo_ids"]),(25,115))
        self.assertEqual((r["gap_one_photo_ids"],r["gap_two_photo_ids"]),(39,76))
        self.assertEqual(r["original_10km_gate"],"HOLD_UNCHANGED")
        self.assertTrue(all(x["quality_review_status"]=="UNREVIEWED" for x in r["photo_candidates"]))
        self.assertTrue(all(x["colour_measurable"] is None for x in r["photo_candidates"]))
        self.assertEqual(len({x["photo_id"] for x in r["photo_candidates"]}),115)

    def test_photo_duplicate_rejected(self):
        b = copy.deepcopy(self.gap1)
        positive = [x for x in b if x["query_status"]=="ELIGIBLE_PUBLIC_METADATA_CANDIDATE_EXISTS"]
        positive[1]["eligible_photo_observation_ids"][0]["photo_id"] = positive[0]["eligible_photo_observation_ids"][0]["photo_id"]
        with self.assertRaisesRegex(ValueError,"Duplicate"):
            assemble(b,self.gap2,self.status,self.queue)

    def test_year_month_mismatch_rejected(self):
        b = copy.deepcopy(self.gap2)
        taxon = next(str(x["inat_taxon_id"]) for x in self.status if x["status"]=="POSSIBLE_COMPLETE_METADATA_ONLY")
        item=next(x for x in b if str(x["inat_taxon_id"])==taxon)
        item["month"]=13
        with self.assertRaisesRegex(ValueError,"month"):
            assemble(self.gap1,b,self.status,self.queue)

    def test_missing_yearly_observer_rejected(self):
        b=copy.deepcopy(self.gap2)
        taxon=next(str(x["inat_taxon_id"]) for x in self.status if x["status"]=="POSSIBLE_COMPLETE_METADATA_ONLY")
        item=next(x for x in b if str(x["inat_taxon_id"])==taxon)
        item["metadata_candidates"]=[]
        with self.assertRaisesRegex(ValueError,"Insufficient"):
            assemble(self.gap1,b,self.status,self.queue)

    def test_no_year_reassignment_or_replacement(self):
        b=copy.deepcopy(self.queue)
        tid=str(next(x for x in self.gap1 if x["query_status"]=="ELIGIBLE_PUBLIC_METADATA_CANDIDATE_EXISTS")["inat_taxon_id"])
        row=next(x for x in b if x["cohort"]=="validation" and int(x["priority_tier"])==1 and str(x["inat_taxon_id"])==tid)
        row["target_year"]="1910"
        with self.assertRaises(ValueError):
            assemble(self.gap1,self.gap2,self.status,b)

if __name__=="__main__":
    unittest.main()
