package cn.edu.thu.imr.util;

import cn.edu.thu.imr.entity.TimePoint;
import cn.edu.thu.imr.entity.TimeSeries;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
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
   * @param filename
   * @param index
   *          which column besides timestamp should be read
   * @return
   */
  public TimeSeries readData(String filename, int index) {
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
  public ArrayList<Boolean> readLabel(String filename, int index) {
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
   * s
   * @param truthSeries
   * @param resultSeries
   * @param labelList
   * @return
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
