from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

for index, path in list(enumerate(sys.path)):
    if "site-packages" in path.lower():
        sys.path.append(sys.path.pop(index))
