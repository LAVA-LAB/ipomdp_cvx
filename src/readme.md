
# Usage
By default, the SCP method will be used.

```
$ python3 main.py --prism [path to prism file] --intervals [path to intervals file] --spec [spec as string] --type ["rew" for reward spec, "prob" for reachability] --threshold [float / int, threshold for termination] 
```
Alternatively, use 
```
--drn [path to drn file]
```
to specificy the model as drn.
Note that either --prism or --drn has to be used, and the other should not be specified.


For instance, to run the satellite model with reward minimization, use:
```
$ python3 main.py --drn "models/satellite/pomdp_attempt_prob_rew_36_sat_065_dist_5_obs_diff_orb_len_1.drn" --intervals "models/satellite/satellite_rew_robust.intervals" --spec "Rmin=?[F \"goal\"]" --type "rew" --threshold 1000 
```

To run the aircraft model with reachability maximization, use
```
$ python3 main.py --drn "models/aircraft/collision_partial_obs_2d_upd_hobs_20_small.drn" --intervals "models/aircraft/collision_partial_obs_2d_upd_hobs_20_small.intervals" --spec "Pmax=?[F \"goal\"]" --threshold 0.75 

```
