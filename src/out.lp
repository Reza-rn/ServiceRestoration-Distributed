\ Model DSR
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
 [ 2 x ^2 + 4 y ^2 + 6 z ^2 ] / 2 
Subject To
 R0: x <= 5
 qc0: [ x ^2 + 2 y ^2 + 3 z ^2 ] = 0
Bounds
 y >= -5
 z >= 2
End
