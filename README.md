# covid_rt_comparison

Report in report.pdf

```bash
.
├── README.md
├── constant
│   └── constant.R
├── data
│   ├── processed
│   │   ├── a.csv
│   │   ├── b.csv
│   │   ├── c.csv
│   │   ├── d.csv
│   │   ├── e.csv
│   │   └── rdate.csv
│   └── raw
│       └── owid_Sep5.csv
├── function
│   ├── add_cyclic.R
│   ├── disc_gamma.R
│   ├── gen_data.R
│   ├── get_iwt.R
│   ├── make_plot.R
│   ├── process_data.R
│   ├── run_trend_filter.R
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
│       ├── a_1day_pois.png
│       ├── a_pois.png
│       ├── b.png
│       ├── b_1day_pois.png
│       └── b_pois.png
├── report.Rmd
├── report.pdf
└── run
    ├── run_EpiEstim.R
    ├── run_EpiInvert.R
    ├── run_EpiNow.R
    └── run_Pl_simple.R

14 directories, 35 files
```
