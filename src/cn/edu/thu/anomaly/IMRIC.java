package cn.edu.thu.anomaly;

import java.util.ArrayList;
import Jama.Matrix;
import cn.edu.thu.anomaly.entity.TimePoint;
import cn.edu.thu.anomaly.entity.TimeSeries;

public class IMRIC extends BaseIMR {

  /**
   * initialize Matrix A and B
   * 
   * @param A
   * @param B
   * @param zs
   */
  private void initMatrix(Matrix A, Matrix B, double[] zs) {
    int size = zs.length;

    // A:p*p, B:p*1
    double[][] a = A.getArray(), b = B.getArray();
    double val = 0;
    int realL = 0; // since the index will be 1 smaller

    // aii
    for (int i = 1; i <= p; ++i) {
      val = 0;
      for (int l = p + 1 - i; l <= size - i; ++l) {
        realL = l - 1;
        val += zs[realL] * zs[realL];
      }
      a[i - 1][i - 1] = val;
    }

    // aij=aji
    for (int i = 1; i <= p; ++i) {
      for (int j = i + 1; j <= p; ++j) {
        val = 0;
        int u = j - i;
        for (int l = p + 1 - i; l <= size - i; ++l) {
          realL = l - 1;
          val += zs[realL] * zs[realL - u];
        }
        a[i - 1][j - 1] = val;
        a[j - 1][i - 1] = val;
      }
    }

    // bi
    for (int i = 1; i <= p; ++i) {
      val = 0;
      int u = i;
      for (int l = p + 1; l <= size; ++l) {
        realL = l - 1;
        val += zs[realL] * zs[realL - u];
      }
      b[i - 1][0] = val;
    }
  }

  /**
   * following formulas (18-20)
   * 
   * @param index
   * @param preVal
   * @param val
   * @param A
   * @param B
   * @param zs
   */
  private void update(int index, double preVal, double val, Matrix A, Matrix B, double[] zs) {
    int size = zs.length; // denoting n+1 since starting with 0
    int n = size - 1;
    int zPos = index + p;
    zs[zPos] = val;

    // A:p*p, B:p*1
    double[][] a = A.getArray(), b = B.getArray();
    double aiiVal = val * val - preVal * preVal;
    double aijVal = val - preVal;
    double addVal = 0;

    // the border will minus 1 since the index
    // n already minus 1
    // aii
    for (int i = 1; i <= p; ++i) {
      if (zPos >= (p + 1 - i) - 1 && zPos <= n - i) {
        a[i - 1][i - 1] += aiiVal;
      }
    }

    // aij
    for (int i = 1; i <= p; ++i) {
      for (int j = i + 1; j <= p; ++j) {
        if (zPos < (p + 1 - j) - 1 || zPos > n - i) {
          addVal = 0;
        } else if (zPos >= (p + 1 - j) - 1 && zPos < (p + 1 - i) - 1) {
          addVal = aijVal * zs[zPos + j - i];
        } else if (zPos > n - j && zPos <= n - i) {
          addVal = aijVal * zs[zPos - j + i];
        } else {
          addVal = aijVal * (zs[zPos - j + i] + zs[zPos + j - i]);
        }
        a[i - 1][j - 1] += addVal;
        a[j - 1][i - 1] += addVal;
      }
    }

    // bi
    for (int i = 1; i <= p; ++i) {
      int u = i;
      if (zPos < (p + 1 - u) - 1) {
        addVal = 0;
      } else if (zPos >= (p + 1 - u) - 1 && zPos < (p + 1 + u) - 1) {
        addVal = aijVal * zs[zPos + u];
      } else if (zPos > n - u) {
        addVal = aijVal * zs[zPos - u];
      } else {
        addVal = aijVal * (zs[zPos - u] + zs[zPos + u]);
      }
      b[i - 1][0] += addVal;
    }
  }

  private void compute(double[][] x, double[][] y, double[] zs) {
    int rowNum = y.length;

    // begin iteration
    int index = -1;
    Matrix xMatrix = new Matrix(x);
    Matrix yMatrix = new Matrix(y);
    int iterationNum = 0;
    double val = 0, preVal = 0;

    Matrix aMatrix = new Matrix(p, p), bMatrix = new Matrix(p, 1);
    initMatrix(aMatrix, bMatrix, zs);

    while (true) {
      iterationNum++;

      // Matrix phi = learnParamsOLS(xMatrix, yMatrix);
      Matrix phi = learnParamsIC(aMatrix, bMatrix);
      Matrix yhatMatrix = combine(phi, xMatrix);

      index = repairAMin(yhatMatrix, yMatrix);
      if (index == -1)
        break;

      preVal = yMatrix.get(index, 0);
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

      int zPos = index + p;
      zs[zPos] = val;
      update(index, preVal, val, aMatrix, bMatrix, zs);

      if (iterationNum > maxNumIterations)
        break;
    }

    System.out.println("Stop after " + iterationNum + " iterations");
  }

  private void incrementalCompute1(double[][] x, double[][] y, double[] zs) {
    int size = zs.length; // denoting n+1 since starting with 0
    int rowNum = size - p;

    // initial alpha and beta
    double alpha = 0, beta = 0;

    for (int i = 1; i < size - 1; ++i) {
      alpha += zs[i] * zs[i];
      beta += zs[i] * zs[i - 1];
    }
    alpha += zs[0] * zs[0];
    beta += zs[size - 1] * zs[size - 2];
    double[][] phiArray = new double[1][1];

    // begin iteration
    int index = -1;
    Matrix xMatrix = new Matrix(x);
    Matrix yMatrix = new Matrix(y);
    int iterationNum = 0;
    double val = 0, preVal = 0;

    while (true) {
      iterationNum++;
      phiArray[0][0] = beta / alpha;
      Matrix phi = new Matrix(phiArray);
      // Matrix phi = learnParamsOLS(xMatrix, yMatrix);
      Matrix yhatMatrix = combine(phi, xMatrix);

      index = repairAMin(yhatMatrix, yMatrix);
      if (index == -1)
        break;

      preVal = yMatrix.get(index, 0);
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

      // update alpha
      int zPos = index + p;
      zs[zPos] = val;
      if (zPos <= size - 2) {
        alpha = alpha - preVal * preVal + val * val;
      }
      // update beta
      double zPosPreVal = zs[zPos - 1];
      double zPosNextVal = zPos < size - 1 ? zs[zPos + 1] : 0;
      beta = beta + (val - preVal) * (zPosPreVal + zPosNextVal);

      if (iterationNum > maxNumIterations)
        break;
    }

    System.out.println("Stop after " + iterationNum + " iterations");
  }

  private void incrementalCompute2(double[][] x, double[][] y, double[] zs) {
    int size = zs.length; // denoting n+1 since starting with 0
    int rowNum = size - p;

    // initial alpha, beta, gamma
    double alpha = 0, beta = 0, gamma = 0;

    // [3,n-2] -> [2,size-3]
    for (int i = 2; i < size - 2; ++i) {
      alpha += zs[i] * zs[i];
      beta += zs[i] * zs[i - 1];
      gamma += zs[i] * zs[i - 2];
    }

    double alpha1 = alpha + zs[0] * zs[0];
    double alpha2 = alpha + zs[size - 2] * zs[size - 2];
    double beta2 = beta + zs[1] * zs[0];
    double beta3 = beta + zs[size - 1] * zs[size - 2];
    double gamma3 = gamma + zs[size - 2] * zs[size - 4] + zs[size - 1] * zs[size - 3];
    double det = alpha2 * alpha1 - beta2 * beta2;

    double alphachange = 0;

    double[][] phiArray = new double[2][1];

    // begin iteration
    int index = -1;
    Matrix xMatrix = new Matrix(x);
    Matrix yMatrix = new Matrix(y);
    int iterationNum = 0;
    double val = 0, preVal = 0;

    while (true) {
      iterationNum++;
      phiArray[0][0] = (beta3 * alpha2 - gamma3 * beta2) / det;
      phiArray[1][0] = (-beta3 * beta2 + gamma3 * alpha1) / det;
      Matrix phi = new Matrix(phiArray);
      // Matrix phi = learnParamsOLS(xMatrix, yMatrix);
      Matrix yhatMatrix = combine(phi, xMatrix);

      index = repairAMin(yhatMatrix, yMatrix);
      if (index == -1)
        break;

      preVal = yMatrix.get(index, 0);
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

      // update alpha
      int zPos = index + p;
      zs[zPos] = val;
      alphachange = -preVal * preVal + val * val;
      if (zPos <= size - 2 && zPos >= 2) {
        alpha2 += alphachange;
      }
      if (zPos <= size - 3 && zPos >= 1) {
        alpha1 += alphachange;
      }

      // update beta
      double zPosPreValb2 = zs[zPos - 1];
      double zPosPreValb3 = zPos >= 2 ? zs[zPos - 1] : 0;
      double zPosNextValb2 = zPos <= size - 2 ? zs[zPos + 1] : 0;
      double zPosNextValb3 = zPos <= size - 2 ? zs[zPos + 1] : 0; // not to
                                                                  // exceed size

      beta2 = beta2 + (val - preVal) * (zPosPreValb2 + zPosNextValb2);
      beta3 = beta3 + (val - preVal) * (zPosPreValb3 + zPosNextValb3);
      // update gamma
      double zPosPreValg3 = zPos >= 2 ? zs[zPos - 2] : 0;
      double zPosNextValg3 = zPos <= size - 3 ? zs[zPos + 2] : 0;

      gamma3 = gamma3 + (val - preVal) * (zPosPreValg3 + zPosNextValg3);

      // update det
      det = alpha2 * alpha1 - beta2 * beta2;

      if (iterationNum > maxNumIterations)
        break;
    }

    System.out.println("Stop after " + iterationNum + " iterations");
  }

  private void incrementalCompute3(double[][] x, double[][] y, double[] zs) {
    int size = zs.length; // denoting n+1 since starting with 0
    int rowNum = size - p;

    // initial alpha, beta, gamma
    double alpha = 0, beta = 0, gamma = 0, zeta = 0;

    // [4,n-3] -> [3,size-4]
    for (int i = 3; i <= size - 4; ++i) {
      alpha += zs[i] * zs[i];
      beta += zs[i] * zs[i - 1];
      gamma += zs[i] * zs[i - 2];
      zeta += zs[i] * zs[i - 3];
    }

    double alpha1 = alpha + zs[2] * zs[2] + zs[1] * zs[1] + zs[0] * zs[0];
    double alpha2 = alpha + zs[2] * zs[2] + zs[1] * zs[1] + zs[size - 3] * zs[size - 3];
    double alpha3 =
        alpha + zs[2] * zs[2] + zs[size - 3] * zs[size - 3] + zs[size - 2] * zs[size - 2];

    double beta2 = beta + zs[2] * zs[1] + zs[1] * zs[0] + zs[size - 3] * zs[size - 4];
    double beta3 = beta + zs[2] * zs[1] + zs[size - 3] * zs[size - 4] + zs[size - 2] * zs[size - 3];
    double beta4 = beta + zs[size - 3] * zs[size - 4] + zs[size - 2] * zs[size - 3]
        + zs[size - 1] * zs[size - 2];

    double gamma3 =
        gamma + zs[2] * zs[0] + zs[size - 3] * zs[size - 5] + zs[size - 2] * zs[size - 4];
    double gamma4 = gamma + zs[size - 3] * zs[size - 5] + zs[size - 2] * zs[size - 4]
        + zs[size - 1] * zs[size - 3];

    double zeta4 = zeta + zs[size - 3] * zs[size - 6] + zs[size - 2] * zs[size - 5]
        + zs[size - 1] * zs[size - 4];

    double A = alpha2 * alpha1 - beta2 * beta2;
    double B = -(beta3 * alpha1 - beta2 * gamma3);
    double C = beta3 * beta2 - alpha2 * gamma3;
    double D = -(beta3 * alpha1 - gamma3 * beta2);
    double E = alpha3 * alpha1 - gamma3 * gamma3;
    double F = -(alpha3 * beta2 - beta3 * gamma3);
    double G = beta3 * beta2 - gamma3 * alpha2;
    double H = -(alpha3 * beta2 - gamma3 * beta3);
    double I = alpha3 * alpha2 - beta3 * beta3;
    double det = alpha3 * A + beta3 * B + gamma3 * C;

    double alphachange = 0;

    double[][] phiArray = new double[3][1];

    // begin iteration
    int index = -1;
    Matrix xMatrix = new Matrix(x);
    Matrix yMatrix = new Matrix(y);
    int iterationNum = 0;
    double val = 0, preVal = 0;

    while (true) {
      iterationNum++;
      phiArray[0][0] = (beta4 * A + gamma4 * D + zeta4 * G) / det;
      phiArray[1][0] = (beta4 * B + gamma4 * E + zeta4 * H) / det;
      phiArray[2][0] = (beta4 * C + gamma4 * F + zeta4 * I) / det;
      Matrix phi = new Matrix(phiArray);
      // Matrix phi = learnParamsOLS(xMatrix, yMatrix);
      Matrix yhatMatrix = combine(phi, xMatrix);

      index = repairAMin(yhatMatrix, yMatrix);
      if (index == -1)
        break;

      preVal = yMatrix.get(index, 0);
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

      // update alpha
      int zPos = index + p;
      zs[zPos] = val;
      alphachange = -preVal * preVal + val * val;
      if (zPos <= size - 2 && zPos >= 3) {
        alpha3 += alphachange;
      }
      if (zPos <= size - 3 && zPos >= 2) {
        alpha2 += alphachange;
      }
      if (zPos <= size - 4 && zPos >= 1) {
        alpha1 += alphachange;
      }

      // update beta
      double zPosPreValb2 = zs[zPos - 1];
      double zPosPreValb3 = zPos >= 2 ? zs[zPos - 1] : 0;
      double zPosPreValb4 = zPos >= 3 ? zs[zPos - 1] : 0;
      double zPosNextValb2 = zPos <= size - 3 ? zs[zPos + 1] : 0;
      double zPosNextValb3 = zPos <= size - 2 ? zs[zPos + 1] : 0;
      // in theory size-1, but zPos+1 will exceed
      double zPosNextValb4 = zPos <= size - 2 ? zs[zPos + 1] : 0;

      beta2 = beta2 + (val - preVal) * (zPosPreValb2 + zPosNextValb2);
      beta3 = beta3 + (val - preVal) * (zPosPreValb3 + zPosNextValb3);
      beta4 = beta4 + (val - preVal) * (zPosPreValb4 + zPosNextValb4);
      // update gamma
      double zPosPreValg3 = zPos >= 2 ? zs[zPos - 2] : 0;
      double zPosPreValg4 = zPos >= 3 ? zs[zPos - 2] : 0;
      double zPosNextValg3 = zPos <= size - 3 ? zs[zPos + 2] : 0;
      double zPosNextValg4 = zPos <= size - 3 ? zs[zPos + 2] : 0;

      gamma3 = gamma3 + (val - preVal) * (zPosPreValg3 + zPosNextValg3);
      gamma4 = gamma4 + (val - preVal) * (zPosPreValg4 + zPosNextValg4);
      // update zeta
      double zPosPreValz4 = zPos >= 3 ? zs[zPos - 3] : 0;
      double zPosNextValz4 = zPos <= size - 4 ? zs[zPos + 3] : 0;

      zeta4 = zeta4 + (val - preVal) * (zPosPreValz4 + zPosNextValz4);
      // update others
      A = alpha2 * alpha1 - beta2 * beta2;
      B = -(beta3 * alpha1 - beta2 * gamma3);
      C = beta3 * beta2 - alpha2 * gamma3;
      D = -(beta3 * alpha1 - gamma3 * beta2);
      E = alpha3 * alpha1 - gamma3 * gamma3;
      F = -(alpha3 * beta2 - beta3 * gamma3);
      G = beta3 * beta2 - gamma3 * alpha2;
      H = -(alpha3 * beta2 - gamma3 * beta3);
      I = alpha3 * alpha2 - beta3 * beta3;
      det = alpha3 * A + beta3 * B + gamma3 * C;

      if (iterationNum > maxNumIterations)
        break;
    }

    System.out.println("Stop after " + iterationNum + " iterations");
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
  public TimeSeries mainIMRIC(TimeSeries dirtySeries, TimeSeries labelSeries,
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

    switch (p) {
      case 1:
        incrementalCompute1(x, y, zs);
        break;
      case 2:
        incrementalCompute2(x, y, zs);
        break;
      case 3:
        incrementalCompute3(x, y, zs);
        break;
      default:
        compute(x, y, zs);
    }

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
        modify = dirtySeries.getTimeseries().get(i).getVal() + y[i - p][0];
      }

      tp = new TimePoint(timestamp, modify);
      resultSeries.addTimePoint(tp);
    }

    return resultSeries;
  }

  public static void main(String[] args) {
    String inputFileName = "ild3k.data";

    IMRIC imric = new IMRIC();

    TimeSeries dirtySeries = imric.readData(inputFileName, 1);
    TimeSeries labelSeries = imric.readData(inputFileName, 2);
    TimeSeries truthSeries = imric.readData(inputFileName, 3);
    ArrayList<Boolean> labelList = imric.readLabel(inputFileName, 4);

    int p = 3;
    double delta = 0.1;
    int maxNumIterations = 100000;

    TimeSeries resultSeries =
        imric.mainIMRIC(dirtySeries, labelSeries, labelList, p, delta, maxNumIterations);

    double rms = imric.calcRMS(truthSeries, resultSeries);

    System.out.println("RMS error is " + rms);
  }
}
