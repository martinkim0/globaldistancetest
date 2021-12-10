# globaldistancetest
### Data reproducibility
Original data files can be found in the following links:
- 7B3Y: https://www.rcsb.org/3d-view/7B3Y
- 7C53: https://www.rcsb.org/structure/7C53
- 7COT: https://www.rcsb.org/structure/7COT
- 7JTL: https://www.rcsb.org/structure/7JTL

### Code reproducibility
PyMOL 2.5 was used in `gdt_ts.py` and generating figures.

[ColabFold version 1.3](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb?authuser=1#scrollTo=kOblAo-xetgx) was used for predictions. The following parameters were used in the model:
```
use_amber     false

use_templates true

msa_mode      MMseqs2

model_type    auto

pair_mode     unpaired+paired 

num_recycles  3
```