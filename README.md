# non-Hermitian-wave-packet

% Author: Weicen Dong <orange233@sjtu.edu.cn>

% Created Date: 2025/1/21

Aim of this code: reproduce figures in arXiv:2501.12163 (Non-Hermitian wave-packet dynamics and its realization within a non-Hermitian chiral cavity)

My code can be run in MATLAB (a b c d in code are m a n b in our manuscript, respectively):

please run 'twoD.m' to get FIG. 1. 2D wave-packet dynamics with and without Berry curvature (I use code from line 1 to line 153 to do calculation, but it takes time. Then I save the data, so please run start from line 154)

please run 'PBCmol.m' to get FIG. 3. Band dispersion

please run 'Berry.m' to get FIG. 4. Chern number and Berry curvature

please run 'OBCmol' to get FIG. 5. Chiral edge states on zigzag edges

please run 'OBC_gaussian.m' to get FIG. 6. Edge wave-packet dynamics

please run 'PBC_gaussian_kxky.m' to get FIG. 8. Bulk wave-packet dynamics (I use code from line 1 to line 304 to do calculation, but it takes time. Then I save the data, so please run start from line 305, please ignore warning)

If you use this code, please cite it as below: 
@misc{dong2025,
      title={Non-Hermitian wave-packet dynamics and its realization within a non-Hermitian chiral cavity}, 
      author={Weicen Dong and Qing-Dong Jiang and Matteo Baggioli},
      year={2025},
      eprint={2501.12163},
      archivePrefix={arXiv},
      primaryClass={cond-mat.mes-hall},
      url={https://arxiv.org/abs/2501.12163}, 
}
(might be published soon)
