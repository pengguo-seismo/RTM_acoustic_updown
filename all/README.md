# wave_decomp_code
Reverse time migration (RTM) with up/down wavefield decomposition.

This is for RTM using acoustic wave equation. Please also note that this repository also contains code that was used for testing purpose during the source code development. 

Please refer to "Guo, P., & McMechan, G. A. (2020). Up/down image separation in elastic reverse time migration. Pure and Applied Geophysics, 177(10), 4811-4828." for more details. 

To compile RTM with up/down wavefield decomposition, simply open a command window and type:

./compile_rtm_ud_p.com


here I provide up/down wavefield decomposition using four different implementations, which include:

ud_flag == 1 "Shen, P., Albertin, U. (2015). Up-down separation using hilbert transformed source for causal imaging condition. In SEG technical program expanded abstracts 2015, Society of Exploration Geophysicists, pp. 4175–4179."

ud_flag == 2 "Zhao, Y., Zhang, H., Yang, J., & Fei, T. (2018). Reducing artifacts of elastic reverse time migration by the deprimary technique. Geophysics, 83(6), S569–S577."

ud_flag == 3 "Gao, K., Chi, B., Huang, L. (2017). Elastic least-squares reverse time migration with implicit wavefield separation. In SEG technical program expanded abstracts 2017, Society of Exploration Geophysicists, pp. 4389–4394"

ud_flag == 4 "Duveneck, E. (2018). Up/down separation of seismic depth images. Geophysics, 83(5), S375–S385"
