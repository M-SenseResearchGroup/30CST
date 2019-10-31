# 30CST

This repository hosts a matlab script for automated analysis of wearable sensor data recorded by BioStamp nPoint (MC10, Inc., Lexington, MA) wearable inertial sensors secured to the chest and thigh during the 30-second chair stand test (30CST). A description of the metrics computed and their association with established clinical measures and fall history in a smaple of persons with multiple sclerosis is reported in the following manuscript currently submitted for review:

Tulipani, LJ, Meyer, B, Larie, D, Solomon, AJ, McGinnis, RS. Metrics extracted from a single wearable sensor during sit-stand transitions relate to mobility impairment and fall risk in people with multiple sclerosis. Submitted to Gait and Posture 10.31.2019.

To run the script, please download the example dataset from () and unzip it into the same folder as the code in this repo. 

Main entrypoint for the method is Script_30CST.m which can be run directly from the command line in Matlab. Output is a table, written to the command window, that reports the number of reps completed and accelerometer-derived measures of task performance.

Please direct questions to Dr. Ryan S. McGinnis - ryan.mcginnis@uvm.edu
