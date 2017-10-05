# vldb17-anomaly
Code release of "Time Series Data Cleaning: From Anomaly Detection to Anomaly Repairing" (VLDB 17)
The description of code filesare listed below:

- `IMR.java`: Algorithm 1 in the paper. Use IMR algorithm to repair time series with certain lables.
- `TimePoint.java`: the class for TimePoint indicating a time point
- `TimeSeries.java`: the class for TimeSeries indicating a time sereis.

Datasets
----------
We use following two datasets in the paper:

- GPS data with real errors, which will be made public later
- [Intel Lab Data](http://db.csail.mit.edu/labdata/labdata.html) with synthetic errors. One of the example dataset is `data/ild3k.data`


The schema of the data file contains five columns, 

- timestamp: the timestamp of the data
- dirty: the observation
- label: the value after label(if no labels, then it will be hold the observation value)
- truth: the truth
- isLabel(boolean): whether this point is labeled or not

Parameters
----------
The input and output of **IMR** algorithm is:

Method

```
mainIMR(p,delta,maxNumIterations,dirtySeries,labelSeries,labelList)
```

Input:

```
int p = 3
double delta = 0.1
int maxNum = 100000
TimeSeries dirtySeries
TimeSeries labelSeries
ArrayList<Boolean> labelList
```

Output

```
Timeseries resultSeries
```

Library
----------
[jama.jar](http://math.nist.gov/javanumerics/jama/) is used to implement the code

Citation
----------
If you use this code for your research, please consider citing:

    @article{DBLP:journals/pvldb/ZhangS0Y17,
      author    = {Aoqian Zhang and
                   Shaoxu Song and
                   Jianmin Wang and
                   Philip S. Yu},
      title     = {Time Series Data Cleaning: From Anomaly Detection to Anomaly Repairing},
      journal   = {{PVLDB}},
      volume    = {10},
      number    = {10},
      pages     = {1046--1057},
      year      = {2017},
    }