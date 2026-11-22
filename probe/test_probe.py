"""Tests for probe.py (no network). MV2 / Munim, Inc. Apache-2.0. Made with AI agents."""
import io
import unittest
import zipfile

import probe


class TestDescribe(unittest.TestCase):
    def test_kinds(self):
        self.assertEqual(probe.kind_of(b"PK\x03\x04rest"), "zip")
        self.assertEqual(probe.kind_of(b"<!DOCTYPE html><html>"), "html")
        self.assertEqual(probe.kind_of(b'<?xml version="1.0"?><VOTABLE/>'), "votable")
        self.assertEqual(probe.kind_of(b"source_id,g_flux\n1,2\n"), "csv")
        self.assertEqual(probe.kind_of(b"Service Unavailable"), "text")

    def test_zip_members_and_header(self):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as z:
            z.writestr("EPOCH_PHOTOMETRY/a.csv", "source_id,band,time\n1,G,2.0\n")
        d = probe.describe(buf.getvalue())
        self.assertEqual(d["kind"], "zip")
        self.assertEqual(d["members"], ["EPOCH_PHOTOMETRY/a.csv"])
        self.assertEqual(d["first_csv_header"], "source_id,band,time")

    def test_cases(self):
        c = probe.cases()
        self.assertEqual(len(c), 8)
        self.assertIn("DATA_STRUCTURE=DATAMODEL_STANDARD", c[0][3])
        self.assertNotIn("DATA_STRUCTURE", c[4][3])
        self.assertIn("RELEASE=Gaia+DR4", c[-1][3])


if __name__ == "__main__":
    unittest.main()
