# GNSS-Augmentation

In dense urban environments, utilizing the Global Navigation Satellite System (GNSS) for localization presents numerous challenges due to:
  - Signal obstructions caused by tall buildings and vehicles
  - Multipath interference from surrounding structures. 
This study addresses these challenges, particularly in the context of data collected from smartphones equipped with GNSS receivers. Investigating outlier detection methods to enhance the reliability of GNSS data in urban settings highlights the differences between each technique to give insight into the advantages and disadvantages.

Our data consists of 132 satellites and 1631 epochs, collected by driving a vehicle equipped with Xiaomi and Samsung smartphones around urban areas, simulating real-world usage scenarios. 
Outlier detection methods were implemented to identify and mitigate localization errors caused by non-line-of-sight and multipath anomalies. These methods included:
  - Global Test: 
    - Chi-square test 
  - Local Test:
    - Parametric Test: 
        - T-test
    - Non-Parametric Test: 
        - Local Outlier Factor (LOF)
        - Isolation Forest (IF)
        - Isolation Forest - Local Outlier Factor (IF-LOF)

Implementing these methods aimed to provide a comprehensive approach to outlier detection, balancing computational efficiency with accuracy. 
