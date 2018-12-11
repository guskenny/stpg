#!/bin/bash

cmd="python draw_plot.py $1/SEED_0.csv $1/SEED_0_BEST.csv $1/TIME_0.csv $1/WALL_TIME_0.csv $1/RESTARTS_0.csv"

if [ $# -gt 1 ]
then
  if [ $2 = "-g" ]
  then
    echo "Plotting $1 (groups)"
    cmd="$cmd $1/GROUPS.csv"
  else
    echo "Plotting $1 (times)"
  fi
else
  echo "Plotting $1 (times)"
fi

$cmd

