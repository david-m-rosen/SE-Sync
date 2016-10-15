function [edges, readPoses, weights] = DaveReadG2oDataset2D(input_file_g2o, doSortEdges)
%
%% INPUTS:
% input_file_g2o: name of the g2o file to read
% doSortEdges: if true, tries to order the edges as (0,1) (1,2) etc, putting
% odometry first and then loop closures
%
%% OUTPUTS:
% (notation note: nrNode is the number of nodes, m is the number of edges in the graph)
% pose: matrix in Real(nrNode x 3) (initial guess, with pose i arranged in row i)
% edges: matrix in Real(m x 11) with each line describing a relative pose measurement.
% A line of edges looks like: [id1 id2 dx dy dth I11 I12 I13 I22 I23 I33]
% and describes the relative pose [dx dy dth] between node id2 and id1 (order is important)
% and with information matrix I.
%
% *IMPORTANT*: in the process of reading, indices will become 1-based (more
% suitable for matlab), while the indices in the original file are 0-based 
%
%% COMMENTS (ii):
% .g2o files format:
% VERTEX SECTION: first nrNodes lines:
% the i-th line looks like "VERTEX_SE2 id_i x_i y_i th_i" (initial guess of
% node i), with indices starting from 0,
% EDGES SECTION: last m lines:
% each line looks like "EDGE_SE2 id_i id_j dx_ij  dy_ij  dth_ij
% Omega(1,1), Omega(1,2), Omega(1,2), Omega(2,2), Omega(2,3), Omega(3,3)",
% where Omega is the information matrix corresponding to the relative pose
% measurements ij (edges can also start by EDGE_SE2_SWITCHABLE)
%
% Luca Carlone 
% Georgia Institute of Technology
% 22/10/2012

if nargin < 2
  doSortEdges = 1;
end

fid_g2o = fopen(input_file_g2o,'r');

%% Create matlab variables
weights = []; % vector of angular measurement variances
unorderedFlag = 0;    % flag that checks if nodes are ordered (0,1,2,..)
unorderedPoseIds = 0; % flag that checks if edges are ordered (0,1)(1,2) etc + loop closures

nodeId = 0;
edgeId = 0;
tline = fgets(fid_g2o);
while ischar(tline)
  
  foundNode = (isempty(strfind(tline, 'VERTEX_SE2 '))==0 || isempty(strfind(tline, 'VERTEX2 '))==0);
  %% VERTEX_SE2 lines
  if foundNode==1
    nodeId = nodeId+1;
    [name, id, x, y, th] = strread(tline, '%s %d %f %f %f');
    if (id+1 ~= nodeId)
      unorderedPoseIds = 1;
    end
    readPosesAndId(nodeId,:) = [id x y th];
  end
  
  foundEdge = (isempty(strfind(tline, 'EDGE_SE2 '))==0 || isempty(strfind(tline, 'EDGE2 '))==0);
  %% EDGE_SE2 lines
  if foundEdge==1
    edgeId = edgeId+1;
    [name, id1, id2, dx, dy, dth, I11, I12, I13, I22, I23, I33] = strread(tline, '%s %d %d %f %f %f %f %f %f %f %f %f');
    if edgeId < nodeId && (id1 ~= edgeId-1 || id2 ~= edgeId) % only on spanning tree
      unorderedFlag = 1;
    end
    edges(edgeId,:) = [id1+1, id2+1, dx, dy, dth, I11, I12, I13, I22, I23, I33];
    weights(edgeId) = 1/I33; % angular measurement variance
  end
  tline = fgets(fid_g2o);
end

%% Check validity of the poses we read
if ~exist('readPosesAndId','var')
  warning('G2o file %s does not contain VERTICES',input_file_g2o)
  readPoses = [];
else
  if unorderedPoseIds == 1
    disp('readG2oDataset2D: incorrect node ids in g2o file')
    [indices, ordering] = sort(readPosesAndId(:,1));
    if(max(indices - [0:length(indices)-1]') > 0)
      error('readG2oDataset2D: incorrect node ids, even after ordering')
    end
    readPoses  = readPosesAndId(ordering,2:4);
  else
    readPoses  = readPosesAndId(:,2:4);
  end
end


fclose(fid_g2o);


