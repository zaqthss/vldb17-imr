package cn.edu.thu.anomaly;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import Jama.Matrix;
import cn.edu.thu.anomaly.entity.TimePoint;
import cn.edu.thu.anomaly.entity.TimeSeries;

public class BaseIMR {
  public static double MINVAL = Double.MAX_VALUE;
  public static double MAXVAL = -Double.MAX_VALUE;
  public static String PATH = "data/";

  protected ArrayList<Boolean> labelList; // whether the point is labeled
  protected int p; // AR(p) model
  protected double delta; // converge
  protected int maxNumIterations; // max iteration number

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
  protected Matrix learnParamsOLS(Matrix xMatrix, Matrix yMatrix) {
    Matrix phi = new Matrix(p, 1);

    Matrix xMatrixT = xMatrix.transpose();

    Matrix middleMatrix = xMatrixT.times(xMatrix);
    phi = middleMatrix.inverse().times(xMatrixT).times(yMatrix);

    return phi;
  }
  
  protected Matrix learnParamsIC(Matrix A, Matrix B) {
    Matrix phi = new Matrix(p, 1);
    
    phi = A.inverse().times(B);
    
    return phi;
  }

  /**
   * use phi to combine
   */
  protected Matrix combine(Matrix phi, Matrix xMatrix) {
    Matrix yhatMatrix = xMatrix.times(phi);

    return yhatMatrix;
  }

  /**
   * Absolute minimum
   */
  protected int repairAMin(Matrix yhatMatrix, Matrix yMatrix) {
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
  protected TimeSeries readData(String filename, int index) {
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
  protected ArrayList<Boolean> readLabel(String filename, int index) {
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
  protected double calcRMS(TimeSeries truthSeries, TimeSeries resultSeries) {
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
}
