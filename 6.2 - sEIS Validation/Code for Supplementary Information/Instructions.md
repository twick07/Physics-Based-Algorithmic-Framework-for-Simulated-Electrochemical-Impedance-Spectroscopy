**Instructions**



**To get the results in Supplementary Information, Section S1.2.1 - Section S1.2.2:**



1\. Generate PyBaMM Reference sEIS Data



-> Open PyBaMM\_EIS\_LCO

-> Set SoC to SoC of interest by changing the SoC value in line 17. e.g. SoC = 100 for 100% SoC case.

-> Run PyBaMM\_EIS\_LCO



2\. Generate MATLAB sEIS Data



-> Open MATLAB\_sEIS\_Simulator\_LCO

-> Set SoC to SoC of interest by changing the SoC value in line 4. e.g. SoC = 100 for 100% SoC case.



3\. Repeat above for different SoC Cases



-> To get the results for different SoC cases, change the SoC values in the respective python and MATLAB scripts. Note possible SoC Values are SoC = \[40 60 80 100].



**To get the results in Supplementary Information, Section S1.2.3:**



1\. Generate PyBaMM Timed sEIS Reference Data



-> Open PyBaMM\_EIS\_LCO\_timed

-> Run PyBaMM\_EIS\_LCO\_timed



2\. Generate MATLAB Timed sEIS Reference Data



-> Open Solver\_Comparison\_LCO

-> Run Solver\_Comparison\_LCO





**To get the results in Supplementary Information, Section S1.2.4 - Section S1.2.5**



1. Generate ECM Fit Comparison Data



-> Open ECM\_Fitting\_Code\_LCO

-> Set SoC to SoC of interest by changing the SoC value in line 14. e.g. SoC = 100 for 100% SoC case.

-> Run ECM\_Fitting\_Code\_LCO



2\. Repeat above for different SoC Cases



-> To get the results for different SoC cases, change the SoC values in the respective python and MATLAB scripts. Note possible SoC Values are SoC = \[40 60 80 100].

