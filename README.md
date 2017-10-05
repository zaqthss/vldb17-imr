# vldb17-anomaly
Code release of "Time Series Data Cleaning: From Anomaly Detection to Anomaly Repairing" (VLDB 17)

- data file schema
    + String inputFileName = "ild3k.data"
    + timestamp,dirty,label,truth,isLabel(boolean)
- input
    + int p = 3
    + double delta = 0.1
    + int maxNum = 100000
    + TimeSeries dirtySeries
    + TimeSeries labelSeries
    + ArrayList<Boolean> labelList
- output
    + Timeseries resultSeries

- jar
    + jama.jar is used to implement the code
    + http://math.nist.gov/javanumerics/jama/