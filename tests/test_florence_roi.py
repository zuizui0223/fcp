from fcp_pipeline.florence_roi import box_area_fraction, choose_flower_box


def test_multiple_boxes_reject_scene_sized_group_box():
    image_size = (768, 1024)
    giant = [0, 34, 766, 493]  # ~45% of image: the failed Raphanus v1 scene box
    central_flower = [330, 300, 430, 420]
    side_flower = [80, 150, 180, 270]
    chosen = choose_flower_box([giant, side_flower, central_flower], image_size)
    assert chosen != giant
    assert chosen == central_flower
    assert box_area_fraction(chosen, image_size) <= 0.20


def test_single_large_closeup_box_is_preserved():
    # The pilot Ipomoea flower legitimately fills most of the image, so a large box
    # is only rejected when Florence also supplied smaller individual proposals.
    image_size = (1024, 768)
    closeup = [55, 53, 823, 766]
    assert box_area_fraction(closeup, image_size) > 0.20
    assert choose_flower_box([closeup], image_size) == closeup


def test_small_central_box_can_win_over_larger_edge_box():
    image_size = (1000, 1000)
    central = [450, 450, 550, 550]
    edge = [0, 0, 250, 250]
    assert choose_flower_box([edge, central], image_size) == central


def test_relaxed_filter_is_used_when_all_multi_boxes_exceed_twenty_percent():
    image_size = (1000, 1000)
    near_center = [250, 250, 800, 800]  # 30.25%
    edge = [0, 0, 600, 600]  # 36%
    chosen = choose_flower_box([edge, near_center], image_size)
    assert chosen == near_center
    assert box_area_fraction(chosen, image_size) <= 0.35


def test_no_valid_boxes_returns_none():
    assert choose_flower_box([], (100, 100)) is None
    assert choose_flower_box([[10, 10, 5, 20]], (100, 100)) is None
