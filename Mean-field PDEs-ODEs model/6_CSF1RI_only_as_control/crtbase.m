% CRTBASE.m - Create base vector 
%
% This function creates a vector containing the base of the loci
% in a chromosome.
%此函数创建一个包含染色体中基因座碱基的向量
%
% Syntax: BaseVec = crtbase(Lind, Base)
%
% Input Parameters:
%
%		Lind	- A scalar or vector containing the lengths
%		    	  of the alleles.  Sum(Lind) is the length of
%		    	  the corresponding chromosome.
%     包含等位基因长度的标量或向量。 Sum(Lind) 是对应染色体的长度。
%
%		Base	- A scalar or vector containing the base of
%		    	  the loci contained in the Alleles.
%     包含等位基因中所含基因座碱基的标量或向量。
%
% Output Parameters:
%
%		BaseVec	- A vector whose elements correspond to the base
%		    	  of the loci of the associated chromosome structure.
%                 一个向量，其元素对应于相关染色体结构基因座的碱基。
%
% Author: Andrew Chipperfield
% Date: 19-Jan-94
%
% Tested under MATLAB v6 by Alex Shenfield (17-Jan-03)

function BaseVec = crtbase(Lind, Base)
%第一次Lind=152,下面ml=1 LenL=1
[ml LenL] = size(Lind) ;%ml=1；LenL=1；
if nargin < 2 
	Base = 2 * ones(LenL,1) ; % default to base 2
end
[mb LenB] = size(Base) ;

% check parameter consistency
%第一次ml=mb=1 描述的是Lind和Base的行数，LenL=LenB=1 描述的是Lind和Base的列数
if ml > 1 | mb > 1
	error( 'Lind or Base is not a vector') ;
elseif (LenL > 1 & LenB > 1 & LenL ~= LenB) | (LenL == 1 & LenB > 1 ) 
	error( 'Vector dimensions must agree' ) ;
elseif LenB == 1 & LenL > 1
	Base = Base * ones(LenL,1) ;
end

BaseVec = [] ;
for i = 1:LenL  %第一次LenL=1
	BaseVec = [BaseVec, Base(i)*ones(Lind(i),1)'];%等于原来的Base*Lind(=264)最后得到一个1*264的矩阵
end