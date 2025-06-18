% REP.m         Replicate a matrix
%
%是一个复制矩阵的函数，在example可以看出来
% This function replicates a matrix in both dimensions. 
%此函数在两个维度上复制矩阵
%
% Syntax:       MatOut = rep(MatIn,REPN);
%
% Input parameters:
%   MatIn    - Input Matrix (before replicating)
%
%   REPN     - Vector of 2 numbers, how many replications in each dimension
%              REPN(1): replicate vertically   垂直复制
%              REPN(2): replicate horizontally  水平复制
%
%              Example:
%
%              MatIn = [1 2 3]
%              REPN = [1 2]: MatOut = [1 2 3 1 2 3]
%              REPN = [2 1]: MatOut = [1 2 3;
%                                      1 2 3]
%              REPN = [3 2]: MatOut = [1 2 3 1 2 3;
%                                      1 2 3 1 2 3;
%                                      1 2 3 1 2 3]
%
% Output parameter:
%   MatOut   - Output Matrix (after replicating)
%
%
% Author:   Carlos Fonseca & Hartmut Pohlheim
% History:  14.02.94        file created
%           22.01.03        tested under MATLAB v6 by Alex Shenfield

function MatOut = rep(MatIn,REPN)
%若代入的MatIn是个10×1的列向量
%这个函数只和带入第一个参数MatIn的行数和列数相关；
%这个函数只和REPN索引中的前两个数相关，第一个数控制MatIn的行数，第二个数控制MatIn的列数；
% Get size of input matrix
   [N_D,N_L] = size(MatIn);

% Calculate
%行和列是对称的
   Ind_D = rem(0:REPN(1)*N_D-1,N_D) + 1;%REPN=[1 32]在前面只用过两次
   Ind_L = rem(0:REPN(2)*N_L-1,N_L) + 1;%rem表示除后的余数
% Create output matrix
   MatOut = MatIn(Ind_D,Ind_L);%变成了全部是MatIn数值大小的向量
end
  
   
   

% End of function