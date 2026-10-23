"""Synthetic URL/bytes tests: no provider photos or historic colour labels read."""
from __future__ import annotations
from io import BytesIO
from PIL import Image
import pytest

from scripts.acquisition import run_real_source_roi_pilot as pilot
from scripts.acquisition.repair_real_source_roi_pilot import FIX_IDS, FIRST_REAL_CSV_SHA


def test_official_provider_uppercase_jpg_is_an_original_photo_not_a_rejection(monkeypatch):
    buff=BytesIO()
    Image.new('RGB',(8,8),(50,20,100)).save(buff,format='JPEG')
    seen=[]
    class FakeResponse:
        def __init__(self,blob):self._bio=BytesIO(blob)
        def __enter__(self):return self._bio
        def __exit__(self,*args):self._bio.close()
    def fake_urlopen(request, timeout):
        seen.append(request.full_url)
        return FakeResponse(buff.getvalue())
    monkeypatch.setattr(pilot,'urlopen',fake_urlopen)
    for pid in sorted(FIX_IDS):
        url=f'https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.JPG'
        result=pilot.original_photo_bytes(url,pid)
        assert result==buff.getvalue()
    assert len(seen)==2 and all(s.endswith('/large.JPG') for s in seen)


def test_provider_photo_url_repair_does_not_widen_host_or_photo_id(monkeypatch):
    def banned(*args,**kwargs):
        raise AssertionError('Invalid source URL reached external image downloader')
    monkeypatch.setattr(pilot,'urlopen',banned)
    for link,pid in (
        ('https://example.org/photos/840402/large.JPG','840402'),
        ('https://inaturalist-open-data.s3.amazonaws.com/photos/840402/large.JPG','3201197'),
        ('https://inaturalist-open-data.s3.amazonaws.com/photos/840402/large.gif','840402'),
        ('https://inaturalist-open-data.s3.amazonaws.com/photos/840402/large.JPG?old=1','840402'),
    ):
        with pytest.raises(ValueError,match='Real photo URL is not'):
            pilot.original_photo_bytes(link,pid)
    assert len(FIRST_REAL_CSV_SHA)==64
