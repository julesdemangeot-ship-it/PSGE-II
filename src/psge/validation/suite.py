"""PSGE-II v1.2 -- CLI entry point for the validation campaign.

Invoked as: python -m psge.validation.suite

This is a thin wrapper around campaign.run_campaign() that handles the CLI
presentation. The actual campaign logic lives in campaign.py.
"""

from psge.validation.campaign import run_campaign, print_campaign_table


def main():
    """Run the campaign and display results."""
    print("PSGE-II v1.2 -- Intrinsic Validation Campaign")
    all_passed, rows = run_campaign()
    print_campaign_table(all_passed, rows)
    return 0 if all_passed else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
