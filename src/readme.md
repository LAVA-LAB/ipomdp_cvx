
# Usage
By default, the SCP method will be used.

```
$ python3 main.py --prism [path to prism file] --intervals [path to intervals file] --spec [spec as string] --type ["rew" for reward spec, "prob" for reachability] --threshold [float / int, threshold for termination] 
```

For instance, to run the satellite model with reward minimization, use:
```
$ python3 main.py --drn "models/satellite/pomdp_attempt_prob_rew_36_sat_065_dist_5_obs_diff_orb_len_1.drn" --intervals "models/satellite/satellite_rew_robust.intervals" --spec "Rmin=?[F \"goal\"]" --type "rew" --threshold 1000 
```
