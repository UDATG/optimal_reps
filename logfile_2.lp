\ Model model1
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x_pos0 + x_pos1 + x_pos2 + x_pos3 + x_pos4 + x_pos5 + x_pos6 + x_neg0
   + x_neg1 + x_neg2 + x_neg3 + x_neg4 + x_neg5 + x_neg6
Subject To
 constraint_1: x_pos0 - x_neg0 = -1
 constraint_2: x_pos1 - x_neg1 = 1
 constraint_3: x_pos2 - x_neg2 - q1 = -1
 constraint_4: x_pos3 - x_neg3 + q1 = 0
 constraint_5: x_pos4 - x_neg4 - q0 = 1
 constraint_6: x_pos5 - x_neg5 + q0 = 0
 constraint_7: x_pos6 - x_neg6 - q0 - q1 = 0
Bounds
End
