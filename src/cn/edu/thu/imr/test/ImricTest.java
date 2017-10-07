package cn.edu.thu.imr.test;

import cn.edu.thu.imr.IMRIC;
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
public class ImricTest {


  public static void main(String[] args) {
    String inputFileName = "ild3k.data";

    IMRIC imric = new IMRIC();
    Assist assist = new Assist();

    TimeSeries dirtySeries = assist.readData(inputFileName, 1);
    TimeSeries labelSeries = assist.readData(inputFileName, 2);
    TimeSeries truthSeries = assist.readData(inputFileName, 3);
    ArrayList<Boolean> labelList = assist.readLabel(inputFileName, 4);

    int p = 3;
    double delta = 0.1;
    int maxNumIterations = 100000;

    TimeSeries resultSeries =
        imric.mainIMRIC(dirtySeries, labelSeries, labelList, p, delta, maxNumIterations);

    double rms = assist.calcRMS(truthSeries, resultSeries, labelList);

    System.out.println("RMS error is " + rms);
  }

}
