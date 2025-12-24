%% Bias-corrected standart and cr-MREPT -- 11.2025 -- v4
% Works for 7 different simulation phantoms
% Select a phantom, select the desired options and run. 
% 3D cylinder phantom is added.
% Scaling is fixed.
% Resolution phantom is added
clear all; close all;

%% Select phantom
% 1: Uniform 
% 2: 2 Anomalies
% 3: 4 Anomalies
% 4: Brain simulation
% 5: 3D Cylinder
%       5.1: Bottom Slice
%       5.2: Middle Slice
%       5.3: Top Slice
% 6: 3D Brain
% 7: Resolution 
opts.phantom = 2;

%% Other options
% Bias correction
% 0: No bias correction
% 1: Apply ias correction
opts.biasCorrection      = 0;      

% Plot options
% 0: None                           - No graphs are plotted
% 1: Reconstruction results only    - Only conductivity results are plotted
% 2: Results only                   - Obtained results are plotted
% 3: All                            - Sub-results are plotted
% 4: Debug                          - Everything is plotted
opts.plots               = 1;

% Noise 
opts.add_noise_to_data   = 0;
opts.noise_level         = 0.01;

% Diffusion Coefficient 
opts.diff_coef           = -0.01;

% Mesh size
opts.meshSize            = 0.0015;

%% Definitions
constants.mu0     = 4e-7*pi;
constants.f0      = 127e6;
constants.w0      = 2*pi*constants.f0;
constants.eps0    = 8.854e-12;
constants.scl     = 1;

%% Load data and set-up
if(opts.phantom == 1)
    load('sigma_uniform.mat');
    load('B1p_complex_water_bc.mat');
    load('ROI_121_0.15_v2.mat');
    
elseif(opts.phantom == 2)
    load('sigma_2anomalies.mat');
    load('B1p_complex_birdcage_v3.mat');
    load('ROI_121_0.15_v2.mat');
    
elseif(opts.phantom == 3)
    load('sigma_4anomalies.mat');
    load('B1p_complex_birdcage_4nomalies_v1.mat')
    load('ROI_121_0.15_v2.mat')
    
elseif(opts.phantom == 4)
    load('sigma_head_transposed.mat');
    load('head4_sim_b1p.mat')
    load('ROI_161_0.15_v2.mat');
    data.mesh.complex.ch1 = (data.mesh.complex.ch1);
    constants.scl = 0.95;

elseif(opts.phantom == 5.1)
    load('cyl3D_sliceBot_v1.mat');
    load('cyl3D_sigma_slice_V1.mat')
    load('ROI_2D_51x51x51_1mm.mat');

elseif(opts.phantom == 5.2)
    load('cyl3D_slice_v1.mat');
    load('cyl3D_sigma_slice_V1.mat')
    load('ROI_2D_51x51x51_1mm.mat');

elseif(opts.phantom == 5.3)
    load('cyl3D_sliceUp_v1.mat');
    load('cyl3D_sigma_slice_V1.mat')
    load('ROI_2D_51x51x51_1mm.mat');

elseif(opts.phantom == 6)
    load('data_slice.mat');
    load('sigma_slice.mat')
    load('ROI_xy.mat');
    sigma_true = (slice.xy);
    data.mesh.complex.ch1 = transpose(data_slice.xy);
    constants.scl = 0.95;

elseif(opts.phantom == 7)
    load('resolution_Sim_v1.mat');
    load('sigma_resolutionSim.mat')
    load('ROI_131_0.1_v1.mat');
    % data.mesh.complex.ch1 = transpose(data.mesh.complex.ch1);
    % sigma_true = (slice.xz);
    % data.mesh.complex.ch1 = transpose(data_slice.xz);
    % ROI.yLin = ROI.zLin;
    % ROI.yMesh = ROI.zMesh;
end

if(opts.plots)
    figure
    imagesc(flip(sigma_true))
    title('True Sigma Distribution')
    colormap jet; clim([0 1]);
end

%% EM Simulation for B1+ magnitude
% Find outer boundary
boundary.B = bwboundaries(sigma_true);
boundary.ind = cell2mat(boundary.B(1,:));
boundary.boundary = zeros(size(sigma_true));
for i = 1:length(boundary.ind)
    boundary.boundary(boundary.ind(i,1),boundary.ind(i,2)) = 1;
end
if(opts.plots > 1)
    figure
    imagesc(boundary.boundary)
    title('Outer Object Boundary')
end

% Set-up EM simulation
EMmodel.model = createpde("electromagnetic","harmonic");
EMmodel.model.FieldType = "magnetic";
EMmodel.SQ = [3,4,-1,-1,1,1,-1,1,1,-1]';
boundary.indHold = boundary.ind; boundary.ind(end,:) = [];
boundary.aa = boundary.ind*ROI.dx - max(ROI.xLin)-ROI.dx;

EMmodel.P1 = [2,length(boundary.ind),boundary.aa(:)']'; 
EMmodel.mask = [2,length(boundary.ind),0.9*boundary.aa(:)']'; 
EMmodel.SQ = [EMmodel.SQ;zeros(length(EMmodel.P1) - length(EMmodel.SQ),1)];

EMmodel.gm = [EMmodel.SQ,EMmodel.P1]; 
EMmodel.g = decsg(EMmodel.gm);
EMmodel.gg = geometryFromEdges(EMmodel.model,EMmodel.g);

EMmodel.EdgeID_all = faceEdges(EMmodel.gg,1);
EMmodel.EdgeID_object = faceEdges(EMmodel.gg,2);
EMmodel.outer_edge = setdiff(EMmodel.EdgeID_all,EMmodel.EdgeID_object); 

if(opts.plots > 2)
    figure
    pdegplot(EMmodel.model,"FaceLabels","on")
end

EMmodel.model.VacuumPermittivity = constants.mu0;
EMmodel.model.VacuumPermeability = constants.eps0;

electromagneticProperties(EMmodel.model,"Face",1, ...
                                    "RelativePermittivity",1, ...
                                    "RelativePermeability",1, ...
                                    "Conductivity",0);
electromagneticProperties(EMmodel.model,"Face",2, ...
                                    "RelativePermittivity",60, ...
                                    "RelativePermeability",1, ...
                                    "Conductivity",0.5);

EMmodel.H       = @(location,state) [-1;-1i]*1; 
EMmodel.BC      =    electromagneticBC(EMmodel.model,"Edge",EMmodel.outer_edge,"MagneticField",EMmodel.H);
EMmodel.mesh    = generateMesh(EMmodel.model,"Hmax",0.005);
EMmodel.result  = solve(EMmodel.model,"Frequency",constants.w0);

EMmodel.Hp = (EMmodel.result.MagneticField.Hx + 1i*EMmodel.result.MagneticField.Hy)/2;
EMmodel.Hm = conj(EMmodel.result.MagneticField.Hx - 1i*EMmodel.result.MagneticField.Hy)/2;

EMmodel.ROI_sim.ind = abs(EMmodel.model.Mesh.Nodes(1,:))<1.1*max(ROI.xLin) & abs(EMmodel.model.Mesh.Nodes(2,:))<1.1*max(ROI.yLin);
EMmodel.ROI_sim.x   = EMmodel.model.Mesh.Nodes(1,EMmodel.ROI_sim.ind);
EMmodel.ROI_sim.y   = EMmodel.model.Mesh.Nodes(2,EMmodel.ROI_sim.ind);
EMmodel.Hp_cut = EMmodel.Hp(EMmodel.ROI_sim.ind);
EMmodel.Hm_cut = EMmodel.Hm(EMmodel.ROI_sim.ind);

results.Hp_cart = reshape(griddata(EMmodel.ROI_sim.x,EMmodel.ROI_sim.y,EMmodel.Hp_cut,ROI.xLin,ROI.yLin),size(sigma_true,1),size(sigma_true,2));
results.Hm_cart = reshape(griddata(EMmodel.ROI_sim.x,EMmodel.ROI_sim.y,EMmodel.Hm_cut,ROI.xLin,ROI.yLin),size(sigma_true,1),size(sigma_true,2));

if(opts.plots == 4)
    figure
    pdegplot(EMmodel.model,"FaceLabels","on");
    hold on;
    pdeplot(EMmodel.model,"XYData",abs(EMmodel.Hp));
    title('B1+ Magnitude - all')
end

if(opts.plots > 1)
    figure
    imagesc(abs(results.Hp_cart))
    title('B1+ Magnitude - center')
end
%% Create masked ROI and mesh
% Uses outer object boundary automatically

[ROI.masked.ps, ROI.masked.es, ROI.masked.ts]=initmesh(decsg(EMmodel.mask),'hmax',opts.meshSize);
ROI.masked.x = ROI.masked.ps(1,:)*constants.scl;
ROI.masked.y = ROI.masked.ps(2,:)*constants.scl;
ROI.masked.z = zeros(size(ROI.masked.x));
ROI.masked.faces = ROI.masked.ts(1:3,:);

ROI.masked.pp = [ROI.masked.x; ROI.masked.y];

ROI.masked.edgeNodes = unique(ROI.masked.es(1:2,:))';
ROI.masked.innerNodes = ones(size(ROI.masked.x));
ROI.masked.innerNodes(ROI.masked.edgeNodes) = 0;
ROI.masked.innerNodes = find(ROI.masked.innerNodes);
ROI.masked.pp = [ROI.masked.x; ROI.masked.y];
ROI.masked.tt = [ROI.masked.faces;ones(1,size(ROI.masked.faces,2))];

%% Create masked ROI and mesh
% Define ROI manually

% % rcn=12;rows3=zeros(1,rcn)';cols3=rows3;
% % radius=0.055;
% % for i=1:rcn
% %     rows3(i)=radius*cos(pi/6*i);
% %     cols3(i)=radius*sin(pi/6*i);
% % end
% % mask    = poly2mask((cols3+0.09)/0.0015+1,(rows3+0.09)/0.0015+1,121,121);
% % p1  = [2; length(rows3); rows3; cols3];
% 
% cols3 =     transpose([0 -0.0337 -0.0604 -0.0698 -0.0442 -0.0153  ...
%             0  0.0153  0.0442  0.0698  0.0604  0.0337 ]);
% rows3 =    -1* transpose([0.0760  0.0760  0.0465  0.0113 -0.0739 -0.0887 ...
%             -0.0866 -0.0887 -0.0739  0.0113  0.0465  0.0760 ]);
% mesh_size=0.0015;
% p1  = [2; length(rows3);  cols3; rows3]; %%% This line here
% 
% [ROI.ps, ROI.es, ROI.ts]=initmesh(decsg(p1),'hmax',opts.meshSize);
% 
% ROI.masked.x = ROI.ps(1,:);
% ROI.masked.y = ROI.ps(2,:);
% ROI.masked.z = zeros(size(ROI.masked.x));
% ROI.masked.faces = ROI.ts(1:3,:);
% 
% ROI.pp = [ROI.masked.x; ROI.masked.y];
% 
% ROI.masked.edgeNodes = unique(ROI.es(1:2,:))';
% ROI.innerNodes = ones(size(ROI.masked.x));
% ROI.innerNodes(ROI.masked.edgeNodes) = 0;
% ROI.masked.innerNodes = find(ROI.innerNodes);
% ROI.masked.pp = [ROI.masked.x; ROI.masked.y];
% ROI.masked.tt = [ROI.masked.faces;ones(1,size(ROI.masked.faces,2))];

%% Start Reconstruction
tic
recon.ang = (angle(data.mesh.complex.ch1));
recon.abs = (abs(results.Hp_cart));
% load('cube2D_v3.mat')
% recon.abs = (abs(data.mesh.complex.ch1));

if(opts.plots == 4)
    figure
    imagesc(recon.ang)
    title('Phase data')
end

% Add noise
if opts.add_noise_to_data==1
    recon.stdnoise = opts.noise_level;
    recon.HpCnoise = recon.stdnoise * (rand(size(recon.ang)) - 0.5);
    recon.ang = recon.ang + recon.HpCnoise;
end

% Gaussian filter % How to pick this automatically?
% recon.h2 = fspecial('gaussian',3,0.5); 
recon.h2 = fspecial('gaussian',5,2);
% recon.h2 = fspecial('gaussian',[7 7],5); 

recon.ang_filtered = filter2(recon.h2,recon.ang);
recon.abs_filtered = filter2(recon.h2,recon.abs);
recon.abs_log      = log(recon.abs_filtered);

if(opts.plots == 4)
    figure
    imagesc(recon.ang_filtered)
    title('Filtered phase data')
end

[recon.dx, recon.dy] = gradient(recon.ang_filtered,ROI.dx,ROI.dy);
recon.lap =  4*del2(recon.ang_filtered,ROI.dx,ROI.dy);

if(opts.plots == 4)
    figure
    imagesc(recon.lap)
    title('Laplacian of phase data')
end


recon.trimesh.phase.lap    = griddata(ROI.xLin,ROI.yLin,recon.lap(:),ROI.masked.x,ROI.masked.y,'linear');
recon.trimesh.PhaseLapFace = mean(recon.trimesh.phase.lap(ROI.masked.faces),1);%

recon.trimesh.phase.phase   = griddata(ROI.xLin,ROI.yLin,recon.ang(:),ROI.masked.x,ROI.masked.y,'linear');
recon.trimesh.phase.dx      = griddata(ROI.xLin,ROI.yLin,recon.dx(:),ROI.masked.x,ROI.masked.y,'linear');
recon.trimesh.phase.dy      = griddata(ROI.xLin,ROI.yLin,recon.dy(:),ROI.masked.x,ROI.masked.y,'linear');
recon.trimesh.PhasedxFace   = mean(recon.trimesh.phase.dx(ROI.masked.faces),1);
recon.trimesh.PhasedyFace   = mean(recon.trimesh.phase.dy(ROI.masked.faces),1);

if(opts.plots == 4)
    figure
    pdesurf(ROI.masked.pp,ROI.masked.tt,recon.trimesh.phase.phase');
    colormap jet;
    title('Phase data in triangular mesh')
    
    figure
    pdesurf(ROI.masked.pp,ROI.masked.tt,recon.trimesh.phase.dx');
    colormap jet;
    title('x-derivative of phase data in triangular mesh')

    figure
    pdesurf(ROI.masked.pp,ROI.masked.tt,recon.trimesh.phase.dy');
    colormap jet;
    title('y-derivative of phase data in triangular mesh')

    figure
    pdesurf(ROI.masked.pp,ROI.masked.tt,recon.trimesh.phase.lap');
    colormap jet;
    title('Laplacian of phase data in triangular mesh')
end

[recon.complexMagnitude.dx, recon.complexMagnitude.dy] = gradient(recon.abs_log,ROI.dx,ROI.dy);
recon.trimesh.Magnitude.dx     = griddata(ROI.xLin,ROI.yLin,recon.complexMagnitude.dx(:),ROI.masked.x,ROI.masked.y,'linear');
recon.trimesh.Magnitude.dy     = griddata(ROI.xLin,ROI.yLin,recon.complexMagnitude.dy(:),ROI.masked.x,ROI.masked.y,'linear');
recon.trimesh.MagnitudedxFace = mean(recon.trimesh.Magnitude.dx(ROI.masked.faces),1);
recon.trimesh.MagnitudedyFace = mean(recon.trimesh.Magnitude.dy(ROI.masked.faces),1);

recon.absFace = mean(recon.abs_filtered(ROI.masked.faces),1);

recon.mag_term = 2 * (recon.trimesh.MagnitudedxFace .* recon.trimesh.PhasedxFace + recon.trimesh.MagnitudedyFace .* recon.trimesh.PhasedyFace);
recon.mag_term_nodes = 2 * (recon.trimesh.Magnitude.dx .* recon.trimesh.phase.dx + recon.trimesh.Magnitude.dy .* recon.trimesh.phase.dy);

if(opts.plots == 4)
    figure
    imagesc(recon.abs_log);
    colormap jet;
    title('Magnitude data in cartesian mesh')
    
    figure
    pdesurf(ROI.masked.pp,ROI.masked.tt,recon.trimesh.Magnitude.dx');
    colormap jet;
    title('x-derivative of magnitude data in triangular mesh')

    figure
    pdesurf(ROI.masked.pp,ROI.masked.tt,recon.trimesh.Magnitude.dy');
    colormap jet;
    title('y-derivative of magnitude data in triangular mesh')
end

if(opts.plots > 2)
    figure
    pdesurf(ROI.masked.pp,ROI.masked.tt,recon.mag_term);colormap jet;
    title('Magnitude term for bias correction')
end

% Standart
if opts.biasCorrection
    recon.cond.std  = (recon.trimesh.phase.lap + recon.mag_term_nodes)/constants.w0/constants.mu0;
else
    recon.cond.std  = recon.trimesh.phase.lap/constants.w0/constants.mu0;
end

% crMREPT
if opts.biasCorrection
[recon.A,recon.Dt2n] = form_full_matrix_trimesh_with_new_diffusion_matrix...
            (ROI.masked.pp,ROI.masked.tt,recon.trimesh.phase.phase,(recon.trimesh.PhasedxFace),(recon.trimesh.PhasedyFace),recon.trimesh.PhaseLapFace + recon.mag_term,opts.diff_coef,1);
else
[recon.A,recon.Dt2n] = form_full_matrix_trimesh_with_new_diffusion_matrix...
            (ROI.masked.pp,ROI.masked.tt,recon.trimesh.phase.phase,(recon.trimesh.PhasedxFace),(recon.trimesh.PhasedyFace),recon.trimesh.PhaseLapFace,opts.diff_coef,1);
end

recon.b         = constants.w0 * constants.mu0 * ones(1,length(ROI.masked.faces));
recon.b         = recon.b.';
recon.sol       = recon.A\recon.b;
recon.cond.cr    = 1./ recon.sol;

% BC
local.b         = constants.w0 * constants.mu0 * ones(1,length(ROI.masked.faces));
local.b         = local.b.';
% nnn=1:length(ROI_local.ind.all); nnn(ROI_local.outernodes)=-1; 
% midnodes = nnn(nnn>0);
Apn1 = recon.A(:,ROI.masked.innerNodes); 

bc_out = 1/0.3;
bpn1 = local.b - recon.A(:,ROI.masked.edgeNodes)*bc_out*ones(length(ROI.masked.edgeNodes),1);
rrr=Apn1\bpn1;
sol=zeros(length(ROI.masked.x),1);
sol(ROI.masked.edgeNodes) = bc_out;
sol(ROI.masked.innerNodes) = rrr;

local.cond.crBC = 1./sol;
toc

%%
if(opts.plots)
    figure
    pdesurf(ROI.masked.pp,ROI.masked.tt,recon.cond.std');colormap jet;
    colormap jet; clim([0 1]);%zlim([0 1]);%ylim([min(ROI.masked.y) max(ROI.masked.y)]);xlim([min(ROI.masked.x) max(ROI.masked.x)]);
    title('Standart MREPT')

    figure
    pdesurf([ROI.masked.x; ROI.masked.y],ROI.masked.tt,recon.cond.cr+0.06);colormap jet;
    colormap jet; clim([0 1]);%zlim([0 1]);%%ylim([min(ROI.masked.y) max(ROI.masked.y)]);xlim([min(ROI.masked.x) max(ROI.masked.x)]);
    title('cr-MREPT - No BC')

    figure
    pdesurf([ROI.masked.x; ROI.masked.y],ROI.masked.tt,local.cond.crBC);colormap jet;
    colormap jet; clim([0 1]);%zlim([0 1]);%%ylim([min(ROI.masked.y) max(ROI.masked.y)]);xlim([min(ROI.masked.x) max(ROI.masked.x)]);
    title('cr-MREPT - BC:0.3')
end