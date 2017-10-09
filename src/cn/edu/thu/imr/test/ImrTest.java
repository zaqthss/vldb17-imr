package cn.edu.thu.imr.test;

import cn.edu.thu.imr.IMR;
import cn.edu.thu.imr.entity.TimeSeries;
import cn.edu.thu.imr.util.Assist;
import java.util.ArrayList;

/**
 * Created by Stoke on 2017/10/7.
 * E-mail address is zaqthss2009@gmail.com
 * Copyright © Stoke. All Rights Reserved.
 *
 * @author Stoke
 */
public class ImrTest {

  public static void main(String[] args) {
    String inputFileName = "ild3k.data";

    Assist assist = new Assist();
    String splitOp = ",";

    TimeSeries dirtySeries = assist.readData(inputFileName, 1, splitOp);
    TimeSeries labelSeries = assist.readData(inputFileName, 2, splitOp);
    TimeSeries truthSeries = assist.readData(inputFileName, 3, splitOp);
    ArrayList<Boolean> labelList = assist.readLabel(inputFileName, 4, splitOp);

    double rmsDirty = assist.calcRMS(truthSeries, dirtySeries, labelList);
    System.out.println("Dirty RMS error is " + rmsDirty);

    int p = 3;
    double delta = 0.1;
    int maxNumIterations = 100000;

    IMR imr = new IMR();
    TimeSeries resultSeries =
        imr.mainIMR(dirtySeries, labelSeries, labelList, p, delta, maxNumIterations);

    double rms = assist.calcRMS(truthSeries, resultSeries, labelList);

    System.out.println("RMS error is " + rms);
  }

}
