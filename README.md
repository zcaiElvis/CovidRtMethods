# covid_rt_comparison

Comparison of Rt on synthetic dataset (EpiFilter, synthetic data 2a and 2b)
Results in <plot>

```bash
.
├── README.md
├── constant
│   └── constant.R
├── data
│   ├── processed
│   │   ├── a.csv
│   │   ├── b.csv
│   │   └── rdate.csv
│   └── raw
│       └── owid_Sep5.csv
├── function
│   ├── disc_gamma.R
│   ├── gen_data.R
│   ├── get_iwt.R
│   ├── make_plot.R
│   ├── process_data.R
│   └── util.R
├── model
│   ├── pls
│   │   ├── penalties_nonsmooth.R
│   │   └── penalties_smooth.R
│   └── vae
├── plot
│   ├── epiestim
│   │   ├── a.png
│   │   └── b.png
│   ├── epiinvert
│   ├── epinow
│   │   ├── a.png
│   │   └── b.png
│   └── pls
│       ├── a.png
│       ├── a_pois.png
│       ├── b.png
│       └── b_pois.png
└── run
    ├── run_EpiEstim.R
    ├── run_EpiInvert.R
    ├── run_EpiNow.R
    └── run_Pl_simple.R

14 directories, 26 files


```
