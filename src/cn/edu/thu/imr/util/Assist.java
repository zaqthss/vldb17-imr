package cn.edu.thu.imr.util;

import cn.edu.thu.imr.entity.TimePoint;
import cn.edu.thu.imr.entity.TimeSeries;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by Stoke on 2017/10/7.
 * E-mail address is zaqthss2009@gmail.com
 * Copyright © Stoke. All Rights Reserved.
 *
 * @author Stoke
 */
public class Assist {
  public static String PATH = "data/";

  /**
   * Basic attributes: timestamp, dirty, label, truth
   *
   * @param filename filename
   * @param index which column besides timestamp should be read
   * @return data in timeseries form
   */
  public TimeSeries readData(String filename, int index, String splitOp) {
    TimeSeries timeSeries = new TimeSeries();

    try {
      FileReader fr = new FileReader(PATH + filename);
      BufferedReader br = new BufferedReader(fr);

      String line;
      long timestamp;
      double value;
      TimePoint tp;

      while ((line = br.readLine()) != null) {
        String[] vals = line.split(splitOp);
        timestamp = Long.parseLong(vals[0]);
        value = Double.parseDouble(vals[index]);

        tp = new TimePoint(timestamp, value);
        timeSeries.addTimePoint(tp);
      }

      br.close();
      fr.close();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    return timeSeries;
  }

  /**
   *
   * @param filename filename
   * @param index the label index
   * @return labelList
   */
  public ArrayList<Boolean> readLabel(String filename, int index, String splitOp) {
    ArrayList<Boolean> labelList = new ArrayList<>();

    try {
      FileReader fr = new FileReader(PATH + filename);
      BufferedReader br = new BufferedReader(fr);

      String line;
      boolean isLabel;

      while ((line = br.readLine()) != null) {
        String[] vals = line.split(splitOp);
        isLabel = Boolean.parseBoolean(vals[index]);

        labelList.add(isLabel);
      }

      br.close();
      fr.close();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    return labelList;
  }

  /**
   * RMS sqrt(|modify - truth|^2 / len)
   * s
   * @param truthSeries truth
   * @param resultSeries after repair
   * @param labelList labelList
   * @return RMS error
   */
  public double calcRMS(TimeSeries truthSeries, TimeSeries resultSeries, ArrayList<Boolean> labelList) {
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
