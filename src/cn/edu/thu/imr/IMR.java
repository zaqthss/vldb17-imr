package cn.edu.thu.imr;

import java.util.ArrayList;
import Jama.Matrix;
import cn.edu.thu.imr.entity.TimePoint;
import cn.edu.thu.imr.entity.TimeSeries;

public class IMR extends BaseIMR {

  /**
   * 
   * @param dirtySeries dirty
   * @param labelSeries label
   * @param labelList labelList
   * @param p order p
   * @param delta delta
   * @param maxNumIterations the maximum number of iterations allowed
   * @return timeseries after repair
   */
  public TimeSeries mainIMR(TimeSeries dirtySeries, TimeSeries labelSeries,
      ArrayList<Boolean> labelList, int p, double delta, int maxNumIterations) {
    setLabelList(labelList);
    setP(p);
    setDelta(delta);
    setMaxNumIterations(maxNumIterations);

    int size = dirtySeries.getLength();
    int rowNum = size - p;

    // form z
    double[] zs = new double[size];
    for (int i = 0; i < size; ++i) {
      zs[i] =
          labelSeries.getTimeseries().get(i).getVal() - dirtySeries.getTimeseries().get(i).getVal();
    }

    // build x,y for params estimation
    double[][] x = new double[rowNum][p];
    double[][] y = new double[rowNum][1];
    for (int i = 0; i < rowNum; ++i) {
      y[i][0] = zs[p + i];
      for (int j = 0; j < p; ++j) {
        x[i][j] = zs[p + i - j - 1];
      }
    }

    // begin iteration
    int index = -1;
    Matrix xMatrix = new Matrix(x);
    Matrix yMatrix = new Matrix(y);
    int iterationNum = 0;
    double val = 0;

    Matrix phi = null;
    while (true) {
      iterationNum++;

      phi = learnParamsOLS(xMatrix, yMatrix);

      Matrix yhatMatrix = combine(phi, xMatrix);

      index = repairAMin(yhatMatrix, yMatrix);
      if (index == -1) { // converge
        break;
      }

      val = yhatMatrix.get(index, 0);
      // update y
      yMatrix.set(index, 0, val);
      // update x
      for (int j = 0; j < p; ++j) {
        int i = index + 1 + j; // p+i-j-1 \Leftrightarrow p+i = index+p
        if (i >= rowNum)
          break;
        if (i < 0)
          continue;

        xMatrix.set(i, j, val);
      }

      if (iterationNum > this.maxNumIterations) {
        break;
      }
    }

    System.out.println("Stop after " + iterationNum + " iterations");

    // form result series
    TimeSeries resultSeries = new TimeSeries();
    long timestamp;
    double modify;
    TimePoint tp;

    for (int i = 0; i < size; ++i) {
      timestamp = labelSeries.getTimeseries().get(i).getTimestamp();
      if (labelList.get(i)) {
        modify = labelSeries.getTimeseries().get(i).getVal();
      } else {
        modify = dirtySeries.getTimeseries().get(i).getVal() + yMatrix.get(i - p, 0);
      }

      tp = new TimePoint(timestamp, modify);
      resultSeries.addTimePoint(tp);
    }

    return resultSeries;
  }
}
