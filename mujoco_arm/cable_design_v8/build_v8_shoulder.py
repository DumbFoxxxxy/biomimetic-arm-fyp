"""Generate V8 shoulder targets and standalone MuJoCo models."""

from __future__ import annotations

import json

from v8_shoulder import (
    CONTACT_PATH,
    EXPORT_ROOT,
    OPEN_CHAIN_PATH,
    TARGETS_PATH,
    build_shoulder_tree,
    generate_targets,
    load_config,
    write_model,
)


def main() -> None:
    config = load_config()
    targets = generate_targets(config)
    TARGETS_PATH.write_text(
        json.dumps(targets, indent=2) + "\n",
        encoding="utf-8",
    )
    write_model(build_shoulder_tree(config, roller_count=0), OPEN_CHAIN_PATH)
    write_model(build_shoulder_tree(config, roller_count=1), CONTACT_PATH)
    EXPORT_ROOT.mkdir(parents=True, exist_ok=True)
    print(f"wrote {TARGETS_PATH}")
    print(f"wrote {OPEN_CHAIN_PATH}")
    print(f"wrote {CONTACT_PATH}")


if __name__ == "__main__":
    main()

