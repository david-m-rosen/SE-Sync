function [] = writeG2oDataset3D(inputFile, measurements, edges_id, poses)
%
% Reads an input file in g2o format:
% VERTEX_SE3:QUAT id(0-based) x y z qx qy qz qw (position & quaternion)
% EDGE_SE3:QUAT id1 id2 dx dy dz dqx dqy dqz dqw 
% I11 I12 I13 I14 I15 I16
%     I22 I23 I24 I25 I26 
%         I33 I34 I35 I36
%             I44 I45 I46
%                 I55 I56
%                     I66
% File is stored in the following matlab structures:
% - poses (stores the VERTEX): 
%   -- poses(i).R is the 3x3 rotation of the i-th node (indices are +1, 1-based)
%   -- poses(i).t is the 3D position of the i-th node
% - edges_id: is an mx2 matrix such that if the k-th row is [i j], it describes a measurement between i and j
% - measurements (stores the EDGES):
%   -- measurements.between(k).R is the k-th relative rotation measurement
%   -- measurements.between(k).t is the k-th relative translation measurement
%   -- measurements.between(k).Info is the 6x6 covariance of k-th measurement
%
% NOTE: input file has zero-based indices, but everything is converted to 1-based in matlab
% 
% Luca Carlone
% Georgia Institute of Technology
% Aug 6, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid_g2o = fopen(inputFile,'w');

nrPoses =  length(poses);
m = size(edges_id,1);

for i = 1:nrPoses % write the n+1 vertices, indexed from 0
  %[VERTEX_SE3:QUAT, id, x, y, z, qx, qy, qz, qw]
  x = poses(i).t(1); y = poses(i).t(2); z = poses(i).t(3); 
  R = poses(i).R;
  q = rot2quat(R,1e-4);
  if norm(q)>1e-3
    q = q/norm(q);
  else
    error('norm close to zero for unit quaternion (1)')
  end
  qw = q(1); qx = q(2); qy = q(3); qz = q(4);
  id = i-1;
  fprintf(fid_g2o,'VERTEX_SE3:QUAT %d %f %f %f %.7f %.7f %.7f %.7f\n',id, x, y, z, qx, qy, qz, qw);
end

for k=1:m % write the m edges
   id1=edges_id(k,1)-1; % write to 0-based indices
   id2=edges_id(k,2)-1; % write to 0-based indices
   dt = measurements.between(k).t;
   dx = dt(1); dy = dt(2); dz = dt(3); 
   dR = measurements.between(k).R;
   dq = rot2quat(dR);
   if norm(dq)>1e-3
      dq = dq/norm(dq);
   else
     error('norm close to zero for unit quaternion (2)')
   end
   dqw = dq(1); dqx = dq(2); dqy = dq(3); dqz = dq(4);
   I = measurements.between(k).Info;
                       
   fprintf(fid_g2o,'EDGE_SE3:QUAT %d %d   %f %f %f   %.7f %.7f %.7f %.7f   %f %f %f %f %f %f   %f %f %f %f %f   %f %f %f %f   %f %f %f   %f %f   %f\n', ...
     id1, id2, dx, dy, dz, dqx, dqy, dqz, dqw, ...
      I(1,1), I(1,2), I(1,3), I(1,4), I(1,5), I(1,6), ...
              I(2,2), I(2,3), I(2,4), I(2,5), I(2,6), ...
                      I(3,3), I(3,4), I(3,5), I(3,6), ...
                              I(4,4), I(4,5), I(4,6), ...
                                      I(5,5), I(5,6), ...
                                              I(6,6));
end

fclose(fid_g2o);
       
            
 

