t = results.Time;
ref = results.Data(964:2492,1);
perform = results.Data(964:2492,2:5);
[m,v,s, e_array] = array_test(ref, perform);
res = ref + e_array;
s = s*0.02;