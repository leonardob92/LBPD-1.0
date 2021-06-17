function nodeData = remove_source_leakage_b(nodeDataOrig, protocol)
%REMOVE_SOURCE_LEAKAGE correct ROI time-courses for source leakage
%
% NODEDATA = REMOVE_SOURCE_LEAKAGE(NODEDATAORIG, PROTOCOL)
%   produces orthogonalised node time-courses NODEDATA from 
%   uncorrected node time-courses NODEDATAORIG
%
%   PROTOCOL is a string to switch between various all-to-all 
%   orthogonalisation methods for source-spread correction. It can be:
%     'none'          - No correction applied. 
%     'symmetric'     - Apply orthogonalisation on the parcel time-courses.
%                       This produces orthonormal parcel time-courses
%                       which are as close as possible to the original
%                       time-courses.
%     'closest'       - Apply orthogonalisation on the parcel time-courses.
%                       Start as for the symmetric method, then converge to
%                       a (not orthonormal) orthogonal matrix which is as
%                       close as possible to the original time-courses. 
%     'householder'   - Orthogonalise using a more numerically stable
%                       alternative to the Gram-Schmidt process, dealing
%                       with ROI time-courses in a random order. 
%
    

%	Copyright 2014 OHBA
%	This program is free software: you can redirstribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%	
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%	
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 235 $
%	$LastChangedDate: 2014-08-07 22:09:15 +0100 (Thu, 07 Aug 2014) $
%	Contact: giles.colclough 'at' magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough and Stephen Smith.

switch protocol
    case 'none'
        % no orthogonalisation applied to parcel time-courses
        
    case 'closest'
        % finds closest orthogonal set of vectors by applying symmetric 
        % orthogonalisation then iterating to find closest orthogonal matrix
        nodeData = find_closest_orthogonal_matrix(nodeDataOrig,@ROInets.closest_orthogonal_matrix,protocol);
        
    case 'symmetric'
        % finds closest orthonormal matrix
        nodeData = find_closest_orthogonal_matrix(nodeDataOrig,@ROInets.symmetric_orthogonalise,protocol);
        
    case 'householder'
        nodeData = find_orthogonal_matrix_by_householder_reflections(nodeDataOrig);
        
    otherwise
        error([mfilename ':UnrecognisedOrthMethod'], ...
              'Unrecognised parcel orthogonalisation protocol. \n');
end%switch
end%remove_source_leakage
% -------------------------------------------------------------------------
function nodeData = find_closest_orthogonal_matrix(nodeDataOrig,orthogFunction,protocol)
%FIND_CLOSEST_ORTHOGONAL_MATRIX wrapper on orthogonalisation functions
    nParcels = ROInets.rows(nodeDataOrig);

	if isa(nodeDataOrig, 'meeg')
        currentMontage      = nodeDataOrig.montage('getmontage');
        if isempty(currentMontage)
            name = sprintf('%s leakage correction',protocol);
        else
            name = sprintf('%s leakage correction - %s',protocol,currentMontage.name);
        end

		[~, ~, ~, W] = orthogFunction(transpose(nodeDataOrig(:,:)));   
		nodeData     = add_montage(nodeDataOrig, W',name,nodeDataOrig.chanlabels);
	else
		nodeData = transpose(orthogFunction(transpose(nodeDataOrig)));
	end%if

end%find_closest_orthogonal_matrix
% -------------------------------------------------------------------------
function nodeData = find_orthogonal_matrix_by_householder_reflections(nodeDataOrig, rankErrorMessage)
%FIND_ORTHOGONAL_MATRIX_BY_HOUSEHOLDER_REFLECTIONS wrapper on householder
% orthogonalisation

nParcels = ROInets.rows(nodeDataOrig);

% isRankDeficient = rank(nodeDataOrig) < nParcels; %original

isRankDeficient = rank(nodeDataOrig,0) < nParcels; %leonardo.. trying to adjust the tolerance to see if the matrix is actually rank deficient or if it only seems to be because of the matlab calculation of tolerance
%the tolerance can be also 1e-16, not necessarily 0
%the following lines seem to proof that the rank of the matrix can change
%accordingly to the linear independence of the columns of the matrix
%because of the default tolerance calculated by defualt.
%playing around with rank of matrices, tollerance and precision of calculations
% B = [1 2 3 4 5; 2 3 8 1 3; 1e-15 2e-15 5e-15 4e-15 5e-15]
% rank(B)
% rank(B,1e-16)
% B2 = [1 2 3 4 5; 2 3 8 1 2; 1e-15 2e-15 3e-15 4e-15 5e-15]
% rank(B2)
% rank(B2,1e-16)
%that to my understanding should proof that we could get a smaller rank
%than the number of independent columns in the matrix even if they are
%linearly independent from each other
if isRankDeficient
    rankErrorMessage = ['The ROI time-course matrix is not full rank. \n',    ...
                        '    This prevents you from using an all-to-all ',    ...
                        'orthogonalisation method. \n',                       ...
                        '    Your data have rank %d, and you are looking ',   ...
                        'at %d ROIs. \n',                                     ...
                        '    You could try reducing the number of ROIs, or ', ...
                        'using an alternative orthogonalisation method. \n'];
    error('ROInets:RankError',rankErrorMessage,rank(nodeDataOrig), nParcels);
end%if

permutation                     = randperm(nParcels);
permutationInverse(permutation) = 1:nParcels;

if isa(nodeDataOrig, 'meeg')
	% bugger the permutation bit for this faff
	[~, ~, ~, W] = ROInets.householder_orthogonalise(nodeDataOrig(:,:).'); 
	nodeData     = add_montage(nodeDataOrig, W', 'householder',nodeDataOrig.chanlabels);
else
	nodeData = ROInets.householder_orthogonalise(nodeDataOrig(permutation,:).').'; 
	nodeData = nodeData(permutationInverse,:);
end

end%find_orthogonal_matrix_by_householder_reflections

