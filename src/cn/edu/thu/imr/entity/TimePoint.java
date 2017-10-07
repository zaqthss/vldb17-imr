package cn.edu.thu.imr.entity;

public class TimePoint {
  private long timestamp;
  private double val;

  public TimePoint(long timestamp, double val) {
    setTimestamp(timestamp);
    setVal(val);
  }

  public long getTimestamp() {
    return timestamp;
  }

  public void setTimestamp(long timestamp) {
    this.timestamp = timestamp;
  }

  public double getVal() {
    return val;
  }

  public void setVal(double val) {
    this.val = val;
  }
}
