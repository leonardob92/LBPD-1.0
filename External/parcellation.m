classdef parcellation 
	% Class-based representation of parcellations
	%
	% This class is designed to facilitate working with parcels and standard masks
	% Key tasks are
	%
	% - Reading in a parcellation and matching it to a standard mask
	% - Reshaping between volume and matrix representations of the parcellation
	% - Removing overlap and binarizing parcellations
	% - Return parcel representations compatible with MEG-ROI-nets
	% - Save nii files with qform/xform headers copied from template files
	% - Display volume figures in fslview or Matlab from parcel data
	%
	% For further help, see the README.md file or osl_example_parcellation.m
	%
	% 
	% KEY PROPERTIES
	%  - p.weight_mask - XYZ x Parcels representation of loaded parcellation
	%  - p.template_mask - The underlying brain image
	%
	% KEY MANIPULATIONS
	% - p.to_vol - Convert a matrix to a 4D volume
	% - p.to_matrix - Convert a 4D volume to a matrix
	% - p.binarize - Return weight mask with no overlap or weights
	% 
	% KEY PLOTTING
	% - p.plot_activation - Display montage in Matlab
	% - p.fslview - Display parcellation in fslview
	% - p.fslview(m) - Convert m to volume and display in fslview
	%
	%
	% March 2017 - Romesh Abeysuriya

	properties
	    weight_mask % XYZ x Parcels matrix of parcel weights
	    labels % Names of each ROI in the parcellation
	end

	% To maintain consistency with the template mask on disk, once the parcellation object
	% is initialized, the template mask can never be changed. This makes manipulating
	% the weight mask in memory much more reliable
	properties(SetAccess = immutable)
		template_mask % Values for each voxel in the template - if this is edited in memory, then the version on disk could be different
		template_coordinates % MNI coordinates for each point in the mask
		template_fname % Filename of standard template, useful for copying qform/sform matrix 
	end

	% These properties are dependent on the weight mask, so are only ever set indirectly
	% They are protected rather than dependent so they are only computed once
	properties(SetAccess = protected)
		is_weighted % 'true' if this is a weighted parcellation
		is_overlapping % 'true' if any voxels belong to more than one parcel
		resolution % Mask resolution (mm)
	end

	properties(Dependent)
		n_parcels % Number of parcels
		n_voxels % Number of voxels
	end

	methods

		function self = parcellation(input_mask,labels,template)
			% Take in a mask or a file storing a mask, convert it to a Parcellation
			%
			% INPUTS
			%	- input_mask specifying the parcels themselves
			%	- labels specifying the ROI names (one for each parcel)
			% 	- template specifying the standard mask (normally automatically inferred - if supplied, it must be the name of a .nii file)
			%	  A nii file is required rather than just a matrix because the xform and other header information is required to be able
			% 	  to write nii files via the parcellation object
			%
			% INPUT_MASK formats that can be loaded
			% - XYZ x Parcels 
			% - Vox x Parcels 
			% - Vox x 1 
			% - XYZ x 1
			% - nii file name (nifti file with one of the above formats)
			% - mat file name (containing variables 'mask' for input mask and 'labels' for labels)
			%
			% See set.weight_mask for how input matrix is reshaped and analyzed
			% Note that the input mask will be clipped to the template (a warning will be shown if this is occurs)
			%
			% LABELS formats that can be loaded
			% - Cell array of strings
			% - mat file name containing variable 'labels'
			% - txt file with each parcel name on a single row
			%
			% If no labels are specified, and the input file is not a mat file, default labels will be generated either by
			% - If a .nii file was loaded, looking for a .txt file with the same file name
			% - Otherwise, placeholder labels 'ROI 1', 'ROI 2' etc. will be generated
			
			if isempty(osldir)
				error('OSL must be initialized because the OSL standard masks are used as templates for parcellations');
			end

			% If file is not present, try looking in the OSL parcellations folder
			if ischar(input_mask) && ~exist(input_mask)
				input_mask = fullfile(osldir,'parcellations',input_mask);
			end

			if nargin < 2 || isempty(labels) 
				labels = [];
			else
				if ~exist(labels)
					labels = fullfile(osldir,'parcellations',labels);
				end

				[~,~,ext] = fileparts(labels);

				switch ext
					case '.txt'
						labels = importdata(labels);
					case '.mat'
						d = load(labels);
						labels = d.labels;
					otherwise
						error(sprintf('Unrecognized label extension ''%s'', must end with .txt or .mat',ext));
				end
			end

			% Load the input file
			if ischar(input_mask) && ~isempty(strfind(input_mask,'nii'))  % nifti input
				
				% See if there is a corresponding label file in the same directory
				% Do this first because input_mask gets converted from filename to mask data
				if isempty(labels)
					if ~isempty(strfind(input_mask,'nii.gz'))
						label_fname = strrep(input_mask,'nii.gz','txt');
					elseif ~isempty(strfind(input_mask,'nii'))
						label_fname = strrep(input_mask,'nii','txt');
					end

					if exist(label_fname,'file')
						labels = importdata(label_fname);
					end
                end
                assert(logical(exist(input_mask,'file')),sprintf('Requested file "%s" could not be found',input_mask))
				input_mask = nii.load(input_mask);

			elseif ischar(input_mask) % matlab input
				d = load(input_mask);
				input_mask = d.mask;

				% If no label file is provided, try and get the labels from the .mat file
				if isempty(labels)
					if isfield(d,'labels')
						labels = d.labels;
					end
				end
			elseif length(input_mask) == 1 % Can enter a spatial resolution to retrieve the whole brain mask i.e. 1 parcel
				[~,~,img] = self.guess_template(input_mask);
				input_mask = +logical(img);
			end

			% Guess and load template if required
			if nargin < 3 || isempty(template) 
				[self.resolution,self.template_fname,self.template_mask] = self.guess_template(input_mask);
			else
				assert(ischar(template) && ~isempty(strfind(template,'nii')),'Template must be a .nii file')
				assert(logical(exist(template,'file')),sprintf('Requested file "%s" could not be found',template))
				self.template_fname = template;
				[self.template_mask,res] = nii.load(template);
				self.resolution = res(1);
			end
			
			self.template_coordinates = osl_mnimask2mnicoords(self.template_fname);

			if ndims(self.template_mask) < 3
				[~,~,mask_1] = self.guess_template(input_mask); % If user gave their own template on a standard grid size
				self.template_mask = matrix2vols(self.template_mask,mask_1);
            end
    
			self.weight_mask = input_mask;

			% Assign default labels if required
			if isempty(labels)
				labels = arrayfun(@(x) sprintf('ROI %d',x),1:self.n_parcels,'UniformOutput',false);
			else
				assert(iscell(labels),'parcellation:invalid_labels','Manually specified labels must be a cell array of ROI names or a file name');
				assert(isvector(labels),'parcellation:invalid_labels','Labels must be a cell vector of strings, not a matrix');
				assert(length(labels)==self.n_parcels,'parcellation:invalid_labels',sprintf('Must have one manually specified label for each parcel (%d provided, %d required)',length(labels),self.n_parcels));
			end

			self.labels = labels(:);
		end

		function self = set.weight_mask(self,mask)
			% Possible mask sizes
			% - Vox x 1 
			% - Vox x Parcels 
			% - XYZ x 1 (3D matrix, unweighted non-overlapping parcel if integers OR single weighted non-overlapping parcel if values are not integers)
			% - XYZ x Parcels (4D matrix, unweighted if all values 0 or 1, overlapping if multiple assingmnet)

			% self.template_coordinates is guaranteed to be set in the constructor
			% and is assigned prior to self.weight_mask. Therefore, if we are here 
			% and self.template_coordinates is empty, this function is being called
			% as part of the default loadobj() i.e. loading from disk. In which case,
			% it should be fine to directly use the mask 
			if isempty(self.template_coordinates)
				self.weight_mask = mask;
				return
			end

			% First check the mask is valid
			% If given a 1D or 2D matrix, make sure volumns correspond to voxels
			if isvector(mask)
				mask = mask(:);
			elseif ndims(mask) == 2 && size(mask,1) ~= self.n_voxels && size(mask,2) == self.n_voxels
				mask = mask.';
			end

			% If given a 1D or 2D matrix, make sure the number of voxels is correct
			if ndims(mask) <= 2
				assert(size(mask,1) == self.n_voxels,'parcellation:invalid_mask',sprintf('Number of voxels in parcellation (%d) does not match the mask (%d)',size(mask,1),self.n_voxels));
			else
				% If given 3D or 4D matrix, make sure that the parcel does not have any voxels that aren't part of the mask
				s = size(mask);
				assert(all(s(1:3) == size(self.template_mask)),'parcellation:invalid_mask',sprintf('Parcellation size (%dx%dx%d) is different to template size (%dx%dx%d)',s(1),s(2),s(3),size(self.template_mask,1),size(self.template_mask,2),size(self.template_mask,3)));
			end

			% Expand all masks to XYZ 
			if ndims(mask) < 3
				mask = matrix2vols(mask,self.template_mask);
			end

			if ndims(mask) == 3
				all_integers = all(mod(mask(:),1)==0); % Do we only have integers?
				vals = unique(mask(:));
				is_weighted = ~(all_integers && all(diff(vals) == 1)); % If we read in 3 dimensions, continuous integers mean interpreted parcel index rather than weight
				
				% Expand XYZ x 1 to XYZ x Parcels if we are confident that this is not a single weighted parcel
				if ~is_weighted
					mask = self.integers_to_masks(mask);
				end
			elseif ndims(mask) == 4
				is_weighted = ~all(ismember(mask(:),[0 1])); % If we read in 4 dimensions, then any non-binary value means its weighted
			end

			% Clip the mask to the template
			to_remove = bsxfun(@and,logical(mask),~logical(self.template_mask)); % Indices of voxels to remove from the parcellation because they are outside the mask
			if any(to_remove(:))
				fprintf(2,'Parcellation is being clipped to template\n');
				mask(to_remove) = 0;
			end

			% Finally, check if overlapping
			parcels_per_voxel = sum(mask~=0,4);
			if max(parcels_per_voxel(:)) > 1
				is_overlapping = true;
			else
				is_overlapping = false;
			end

			self.weight_mask = mask;
			self.is_weighted = is_weighted;
			self.is_overlapping = is_overlapping;
		end

	    function parcelflag = parcelflag(self,binarize)
			% Return the parcelflags for ROInets.get_node_tcs
			if nargin < 2 || isempty(binarize) 
				binarize = false;
			end
			if binarize
				p = self.binarize;
			else
				p = self.weight_mask;
			end

			parcelflag = self.to_matrix(p);
		end

		function dat4 = to_vol(self,dat2);
			% Convert Vox x Frames or Parcels x Frames, to XYZ x Frames
			%
			% Commonly used to produce a matrix suitable for use with fslview

			% If the user passed in an XYZ matrix already, check it is a valid volume and return, otherwise throw an error
			s = size(dat2);
			if ndims(s) > 2 
				assert(all(s(1:3) == size(self.template_mask)),'parcellation:to_vol:invalid_size',sprintf('Passed in a matrix with size (%dx%dx%d) but this does not match the template size (%dx%dx%d)',s(1),s(2),s(3),size(self.template_mask,1),size(self.template_mask,2),size(self.template_mask,3)));
				dat4 = dat2;
				return
			end

			% Priority is Vox x Frames or Parcels x Frames, but if we get
			% Frames x Vox then transpose if Frames ~= n_parcels and if
			% we get Frames x Parcels then transpose if Frames ~= n_voxels
			% The first if statement already catches the exceptions above
			if s(1) ~= self.n_voxels && s(1) ~= self.n_parcels
				if s(2) == self.n_voxels || s(2) == self.n_parcels
					dat2 = dat2.';
				else
					error('parcellation:to_vol:invalid_size',sprintf('Unsupported size (%dx%d) - one dimension must have size %d (number of voxels) or %d (number of parcels).',s(1),s(2),self.n_voxels,self.n_parcels));
				end
			end

			% If we have parcels as the first dimension, we need to expand it to voxels
			% (only works if binarized as otherwise don't assume how to expand with weights)
			if size(dat2,1) == self.n_parcels
				d2 = zeros(self.n_voxels,size(dat2,2));
				m = self.value_vector; % Map voxels to parcels
				for k = 1:size(dat2,2)
					for j = 1:self.n_parcels
						d2(m==j,k) = dat2(j,k);
					end
				end
				dat2 = d2;
			end

			dat4 = matrix2vols(dat2,self.template_mask);
		end

		function dat2 = to_matrix(self,dat4)
			% Convert XYZ x Frames to Vox x Frames
			% If no matrix is provided, it will return the matrix
			% representation of the parcellation
			% Also supports converting Parcels x Frames to Vox x Frames

			if nargin < 2 || isempty(dat4) 
				dat4 = self.weight_mask;
			end

			if ndims(dat4) == 2 && size(dat4,1) == self.n_voxels && size(dat4,2) == self.n_parcels
				dat2 = dat4;
				return
			end

			% A 2D matrix input instead of 4D means we took in Parcels x Frames
			% Use to_vol to convert it to a volume first
			if ndims(dat4) == 2
				dat4 = self.to_vol(dat4);
			end
			
			dat2 = vols2matrix(dat4,self.template_mask);
		end

		function p = value_vector(self)
			% Return a Vox x 1 vector where value is parcel index (binarizes if necessary)
			if self.is_weighted || self.is_overlapping
				fprintf(2,'Warning - parcellation is being binarized\n')
			end
			p = self.to_matrix(self.binarize);
            p = bsxfun(@times,p,1:self.n_parcels);
			p = sum(p,2);
		end

		function weight_mask = remove_weights(self)
			% Return a weight mask with weightings removed, but overlap still permitted
			weight_mask = +(self.weight_mask ~= 0);
		end
		
		function weight_mask = remove_overlap(self)
			% Return a weight mask with no overlap, but original weights preserved
			% e.g. If voxel 1 belongs to parcel 1 (0.5) and parcel 2 (0.25) and it's assigned
			% to parcel 1, then the new weights will be parcel 1 (0.5) and parcel 2 (0)
			weight_mask = self.weight_mask.*self.binarize();
		end

		function weight_mask = binarize(self)
			% Return a weight mask with no weights and no overlap
			% binarize() is equivalent to remove_overlap() followed by remove_weights()
			% WARNING - it is NOT equivalent to remove_weights() followed by remove_overlap()
			p = self.to_matrix(self.weight_mask);
			[~,assignments] = max(p,[],2);
			assignments(all(p==0,2)) = 0; % Voxels x 1 with value indicating parcel
			assignment_vol = self.to_vol(assignments); % XYZ x 1 with values indicating parcel
			weight_mask = self.integers_to_masks(assignment_vol,self.n_parcels); % XYZ x n_parcels
		end

		% PROPERTIES
		function n = get.n_parcels(self)
			n = size(self.weight_mask,4);
		end

		function n = get.n_voxels(self)
			n = size(self.template_coordinates,1);
		end

		function plot(self)
			% Standard wrapper for graphical display - calls show_parcellation()
			self.show_parcellation();
		end

		function fslview(self,activation,clim)
			% Display parcellation or parcel activity using fslview
			% 
			% INPUTS
			% - activation - Activation input in volume format or format supported by p.to_vol()
			% 				 By default, the weight mask is used
			% - clim - Colour axis limit, uses fslview default if not specified (auto scaled)
			% 
			% Other examples
			%	p.fslview() - Shows each parcel as a separate volume
			%   p.fslview(p.value_vector) - Shows parcellation as a single volume

			if nargin < 2 || isempty(activation) 
				activation = self.weight_mask;
			end

			if ndims(activation) < 3
				activation = self.to_vol(activation);
			end

			if nargin < 3 || isempty(clim) 
				clim = [min(activation(:)) max(activation(:))];
			end
			
			fname = self.savenii(activation);
			fslview(fname,clim,[],self.template_fname);
			pause(5);
			delete(fname)
		end

		function fsleyes(self,activation,clim)
			% Display parcellation or parcel activity using fsleyes
			% 
			% INPUTS
			% - activation - Activation input in volume format or format supported by p.to_vol()
			% 				 By default, the weight mask is used
			% - clim - Colour axis limit, uses fsleyes default if not specified (auto scaled)
			% 
			% Other examples
			%	p.fsleyes() - Shows each parcel as a separate volume
			%   p.fsleyes(p.value_vector) - Shows parcellation as a single volume

			if nargin < 2 || isempty(activation) 
				activation = self.weight_mask;
			end

			if ndims(activation) < 3
				activation = self.to_vol(activation);
			end

			if nargin < 3 || isempty(clim) 
				clim = [min(activation(:)) max(activation(:))];
			end
			
			fname = self.savenii(activation);
			fsleyes(fname,clim,[],self.template_fname);
			pause(5);
			delete(fname)
		end

		function o = osleyes(self,activation,options)
			if nargin < 3 || isempty(options) 
				options = struct;
			end
			
			if nargin < 2 || isempty(activation) 
				activation = self.weight_mask;
			end
			
			fname = self.savenii(activation);
			o = osleyes({self.template_fname,fname},options);
			delete(fname)
		end


		function output_fname = savenii(self,data,fname)
			% Save a nii file, with qform/xform copied from the original mask file
			if nargin < 3 || isempty(fname) 
				fname = [tempname('.') '.nii.gz'];
			end

			if nargin < 2 || isempty(data) 
				data = self.weight_mask;
			end

			if ischar(data)
				error('parcellation:savenii:invalid_data','Input data was a string. Correct usage is ''savenii(data,fname)''')
			end
			
			if ndims(data) < 3
				data = self.to_vol(data);
			end

			[~,res,xform] = nii.load(self.template_fname);
			output_fname = nii.save(data,res,xform,fname);

		end

		function [coords,weights] = roi_coordinates(self)
			% Retrive the voxel coordinates associated with each ROI
			%
			% Returns a cell array where the cell index is the ROI
			% and the contents is a list of voxel coordinates associated
			% with that ROI

			w = self.to_matrix(self.weight_mask);
			for j = 1:self.n_parcels
				idx = w(:,j) ~= 0;
				coords{j} = self.template_coordinates(idx,:);
				weights{j} = w(idx,j);
			end

		end

		function centers = roi_centers(self)
			% Return the centre-of-gravity of each parcel
			%
			% This is obtained by taking the weighted sum over the MNI coordinates of
			% all of the voxels belonging to the parcel
			[coords,weights] = self.roi_coordinates;
			centers = nan(self.n_parcels,3);
			for j = 1:self.n_parcels
				centers(j,:) = sum(bsxfun(@times,coords{j},weights{j}))./sum(weights{j});
			end
		end

	end

	methods (Static)
		function d4 = integers_to_masks(d3,n_parcels)
			% Convert a 3D parcellation where value indicates parcel number
			% to a 4D parcellation
			% 
			% INPUTS
			% - d3 = XYZ x 1 where all values are integers
			% OUTPUTS
			% - d4 = XYZ x n_parcels (binary parcellation)
			
			assert(ndims(d3) == 3,'Input must be XYZ x 1')
			assert(all(mod(d3(:),1)==0),'Input must only contain integers'); 
			
            % If number of parcels is not specified, infer from data
            if nargin < 2 || isempty(n_parcels)
                n_parcels = max(d3(:));
            end
            
			d4 = zeros([size(d3) n_parcels]);
			for j = 1:n_parcels
				d4(:,:,:,j) = d3 == j;
			end
		end

		[spatial_res,mask_fname,mask] = guess_template(m);
		
	end

end
