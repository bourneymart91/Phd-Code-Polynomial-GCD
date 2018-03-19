function Sk = BuildDTQ_log(fx, gx, k)



m = GetDegree(fx);
n = GetDegree(gx);

C1 = BuildDTQ_log_partition(fx, n-k);
C2 = BuildDTQ_log_partition(gx, m-k);

Sk = [C1 C2];

end

