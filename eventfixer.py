import numpy as np

def fixer(q, time, flux, cluster): # future functionality to cut different events for each cluster?

   events = [160.3, 181.6, 214.63, 230.4, 254.5, 291, 322, 352.35, 372.9, 383.8, 396.2, 443.3, 476, 504, 538, 567, 598.5, 630, 656.3, 660.5, 690.6, 718.04, 735, 762, 805, 845, 906.5, 937, 969.3, 1063.6, 1126.5, 1154.5, 1182.5, 1215.5, 1273.5, 1305.5, 1336.5, 1415.5, 1435.8, 1471.6]
   qs = [1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 9, 10, 10, 10, 11, 12, 13, 14, 14, 15, 15, 15, 16, 16, 17] # quarter corresponding to event
   cut = [13.7, 6.4, 2.17, 3.1, 9.5, 2.41, 5, 2.46, 1.3, 2.2, 0.7, 2.2, 2, 2, 12, 3, 3.5, 5, 7.1, 1.14, 2.4, 6.96, 5, 4, 6, 3, 5.5, 3.6, 3, 1.9, 2.4, 3, 3.5, 3.25, 1.5, 3, 3.5, 4.5, 1.2, 8.4] # how much to cut from start of event
   
   # clip individual safe modes
   for i in range(len(qs)):
      if q == qs[i]:
         flux = np.concatenate((flux[np.where(time <= events[i])], flux[np.where(time >= events[i]+cut[i])]))
         time = np.concatenate((time[np.where(time <= events[i])], time[np.where(time >= events[i]+cut[i])]))
      else:
         pass

   return time, flux
