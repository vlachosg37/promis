from pathlib import Path

from promis.workflow import get_workflow_path, resolve_resource_path


def test_default_packaged_resources_resolve_under_workflow_dir(tmp_path) -> None:
    workflow_dir = Path(get_workflow_path()).resolve()

    defaults = [
        "database/MSI_loci_hg38_coordinates_metadata_exonic_chr_rem_artefacts.csv",
        "database/cytoBand_hg38.txt",
        "scripts",
    ]

    for default_relative in defaults:
        resolved = resolve_resource_path(
            default_relative,
            default_relative,
            str(workflow_dir / default_relative),
            run_dir=tmp_path,
        )

        assert Path(resolved) == (workflow_dir / default_relative).resolve()


def test_custom_relative_resources_resolve_from_run_dir(tmp_path) -> None:
    workflow_dir = Path(get_workflow_path()).resolve()

    cases = [
        (
            "custom/test_loci.csv",
            "database/MSI_loci_hg38_coordinates_metadata_exonic_chr_rem_artefacts.csv",
        ),
        ("custom/cytoBand.txt", "database/cytoBand_hg38.txt"),
        ("custom/scripts", "scripts"),
    ]

    for custom_relative, default_relative in cases:
        resolved = resolve_resource_path(
            custom_relative,
            default_relative,
            str(workflow_dir / default_relative),
            run_dir=tmp_path,
        )

        assert Path(resolved) == (tmp_path / custom_relative).resolve()


def test_absolute_resources_remain_unchanged(tmp_path) -> None:
    workflow_dir = Path(get_workflow_path()).resolve()
    absolute = tmp_path / "custom" / "test_loci.csv"

    resolved = resolve_resource_path(
        str(absolute),
        "database/MSI_loci_hg38_coordinates_metadata_exonic_chr_rem_artefacts.csv",
        str(
            workflow_dir
            / "database/MSI_loci_hg38_coordinates_metadata_exonic_chr_rem_artefacts.csv"
        ),
        run_dir=tmp_path / "other",
    )

    assert Path(resolved) == absolute
