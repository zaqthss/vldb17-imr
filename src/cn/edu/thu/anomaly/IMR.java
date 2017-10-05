package cn.edu.thu.anomaly;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import Jama.Matrix;
import cn.edu.thu.anomaly.entity.TimePoint;
import cn.edu.thu.anomaly.entity.TimeSeries;

public class IMR {
  public static double MINVAL = Double.MAX_VALUE;
  public static double MAXVAL = -Double.MAX_VALUE;
  public static String PATH = "data/";

  private ArrayList<Boolean> labelList; // whether the point is labeled
  private int p; // AR(p) model
  private double delta; // converge
  private int maxNumIterations; // max iteration number

  public void setLabelList(ArrayList<Boolean> labelList) {
    this.labelList = labelList;
  }

  public void setP(int p) {
    this.p = p;
  }

  public void setDelta(double delta) {
    this.delta = delta;
  }

  public void setMaxNumIterations(int maxNumIterations) {
    this.maxNumIterations = maxNumIterations;
  }

  /**
   * learn the parameter can use incremental computation for x and y
   */
  private Matrix learnParamsOLS(Matrix xMatrix, Matrix yMatrix) {
    Matrix phi = new Matrix(p, 1);

    Matrix xMatrixT = xMatrix.transpose();

    Matrix middleMatrix = xMatrixT.times(xMatrix);
    phi = middleMatrix.inverse().times(xMatrixT).times(yMatrix);

    return phi;
  }

  /**
   * use phi to combine
   */
  private Matrix combine(Matrix phi, Matrix xMatrix) {
    Matrix yhatMatrix = xMatrix.times(phi);

    return yhatMatrix;
  }

  /**
   * Absolute minimum
   */
  private int repairAMin(Matrix yhatMatrix, Matrix yMatrix) {
    int rowNum = yhatMatrix.getRowDimension();
    Matrix residualMatrix = yhatMatrix.minus(yMatrix);

    double aMin = MINVAL;
    int targetIndex = -1;
    double yhat, yhatabs;

    for (int i = 0; i < rowNum; ++i) {
      if (labelList.get(i + p)) {
        continue;
      }
      if (Math.abs(residualMatrix.get(i, 0)) < delta) {
        continue;
      }

      yhat = yhatMatrix.get(i, 0);
      yhatabs = Math.abs(yhat);

      if (yhatabs < aMin) { // no need to > 0
        aMin = yhatabs;
        targetIndex = i;
      }
    }

    return targetIndex;
  }

  /**
   * Basic attributes: timestamp, dirty, label, truth
   * 
   * @param filename
   * @param index
   *          which column besides timestamp should be read
   * @return
   */
  private TimeSeries readData(String filename, int index) {
    TimeSeries timeSeries = new TimeSeries();

    try {
      FileReader fr = new FileReader(PATH + filename);
      BufferedReader br = new BufferedReader(fr);

      String line = null;
      long timestamp;
      double value;
      TimePoint tp = null;

      while ((line = br.readLine()) != null) {
        String[] vals = line.split(",");
        timestamp = Long.parseLong(vals[0]);
        value = Double.parseDouble(vals[index]);

        tp = new TimePoint(timestamp, value);
        timeSeries.addTimePoint(tp);
      }

      br.close();
      fr.close();
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    return timeSeries;
  }

  /**
   * 
   * @param filename
   * @param index
   * @return
   */
  private ArrayList<Boolean> readLabel(String filename, int index) {
    ArrayList<Boolean> labelList = new ArrayList<>();

    try {
      FileReader fr = new FileReader(PATH + filename);
      BufferedReader br = new BufferedReader(fr);

      String line = null;
      boolean isLabel;

      while ((line = br.readLine()) != null) {
        String[] vals = line.split(",");
        isLabel = Boolean.parseBoolean(vals[index]);

        labelList.add(isLabel);
      }

      br.close();
      fr.close();
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    return labelList;
  }

  /**
   * RMS sqrt(|modify - truth|^2 / len)
   * 
   * @param timeseries
   * @return
   */
  private double calcRMS(TimeSeries truthSeries, TimeSeries resultSeries) {
    double cost = 0;
    double delta;
    int len = truthSeries.getLength();
    int labelNum = 0;

    for (int i = 0; i < len; ++i) {
      if (labelList.get(i)) {
        labelNum++;
        continue;
      }
      delta = resultSeries.getTimeseries().get(i).getVal()
          - truthSeries.getTimeseries().get(i).getVal();

      cost += delta * delta;
    }
    cost /= (len - labelNum);

    return Math.sqrt(cost);
  }

  /**
   * 
   * @param dirtySeries
   * @param labelSeries
   * @param labelList
   * @param p
   * @param delta
   * @param maxNumIterations
   * @return
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
      zs[i] = labelSeries.getTimeseries().get(i).getVal()
          - dirtySeries.getTimeseries().get(i).getVal();
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
        modify = dirtySeries.getTimeseries().get(i).getVal()
            + yMatrix.get(i - p, 0);
      }

      tp = new TimePoint(timestamp, modify);
      resultSeries.addTimePoint(tp);
    }

    return resultSeries;
  }

  public static void main(String[] args) {
    String inputFileName = "ild3k.data";

    IMR imr = new IMR();

    TimeSeries dirtySeries = imr.readData(inputFileName, 1);
    TimeSeries labelSeries = imr.readData(inputFileName, 2);
    TimeSeries truthSeries = imr.readData(inputFileName, 3);
    ArrayList<Boolean> labelList = imr.readLabel(inputFileName, 4);

    int p = 3;
    double delta = 0.1;
    int maxNumIterations = 100000;

    TimeSeries resultSeries = imr.mainIMR(dirtySeries, labelSeries, labelList,
        p, delta, maxNumIterations);

    double rms = imr.calcRMS(truthSeries, resultSeries);

    System.out.println("RMS error is " + rms);
  }
}
