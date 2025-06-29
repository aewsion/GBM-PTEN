% CRTBP.m - Create an initial population
%
% This function creates a binary population of given size and structure.
%
% Syntax: [Chrom Lind BaseV] = crtbp(Nind, Lind, Base)
%
% Input Parameters:
%
%		Nind	- Either a scalar containing the number of individuals
%	    		  in the new population or a row vector of length two
%			      containing the number of individuals and their length.
%
%		Lind	- A scalar containing the length of the individual
%   			  chromosomes.
%
%		Base	- A scalar containing the base of the chromosome 
%			      elements or a row vector containing the base(s) 
%   			  of the loci of the chromosomes.
%
% Output Parameters:
%
%		Chrom	- A matrix containing the random valued chromosomes 
%			      row wise.
%
%		Lind	- A scalar containing the length of the chromosome.
%
%		BaseV	- A row vector containing the base of the 
%   			  chromosome loci.
%
% Author: Andrew Chipperfield
% Date:	19-Jan-94
%
% Tested under MATLAB v6 by Alex Shenfield (20-Jan-03)

function [Chrom, Lind, BaseV] = crtbp(Nind, Lind, Base)
nargs = nargin ;%送进来多少个参数，那么nargin数值就为多少
%第一次送进来 NIND = 20；Lind=NVAR*PRECI=33*8，Base没有数值

% Check parameter consistency

if nargs >= 1, [mN, nN] = size(Nind) ; end
if nargs >= 2, [mL, nL] = size(Lind) ; end
if nargs == 3, [mB, nB] = size(Base) ; end
%第一次nN=1,不执行nN==2
if nN == 2
   if (nargs == 1) 
      Lind = Nind(2) ; Nind = Nind(1) ; BaseV = crtbase(Lind) ;
   elseif (nargs == 2 & nL == 1) 
      BaseV = crtbase(Nind(2),Lind) ; Lind = Nind(2) ; Nind = Nind(1) ; 
   elseif (nargs == 2 & nL > 1) 
      if Lind ~= length(Lind), error('Lind and Base disagree'); end
      BaseV = Lind ; Lind = Nind(2) ; Nind = Nind(1) ; 
   end
elseif nN == 1
   if nargs == 2
      if nL == 1, BaseV = crtbase(Lind) ;%第一次进入这，得到BaseV是一个1*264，每个数字都是2的向量
      else, BaseV = Lind ; Lind = nL ; end
   elseif nargs == 3
      if nB == 1, BaseV = crtbase(Lind,Base) ; 
      elseif nB ~= Lind, error('Lind and Base disagree') ; 
      else BaseV = Base ; end
   end
else
   error('Input parameters inconsistent') ;
end

% Create a structure of random chromosomes in row wise order, dimensions
% Nind by Lind. The base of each chromosomes loci is given by the value
% of the corresponding element of the row vector base.

Chrom = floor(rand(Nind,Lind).*BaseV(ones(Nind,1),:)) ;
%%这样设置出来的染色体每个位置上的元素全为0或1，因为构造出来的BaseV是一个全是2的矩阵
%X = rand 返回一个在区间 (0,1) 内均匀分布的随机数。
asd=1;%为了断点
end
%Y = floor(X) 将 X 的每个元素四舍五入到小于或等于该元素的最接近整数。
%BaseV(ones(Nind,1),:)是相当于本来BaseV是一个1*264的向量，现在变成了，Nind*264的矩阵了，复制了Nind行
%%每一行和原来的行一样
%最后生成一条染色体，有20*264的0-1的随机数点乘以20*264


% End of file 