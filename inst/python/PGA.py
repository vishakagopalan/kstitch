#from tqdm import tqdm
import numpy as np
import os
#import pyarrow
import pyarrow.parquet as pq
#import matplotlib
#import pyarrow as pa

#import matplotlib.style as mplstyle
#mplstyle.use('fast')

# matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import pandas as pd
#import os
#import scanpy as sc
from multiprocess import Pool
import h5py
from itertools import repeat
# import visvalingamwyatt as vw
#import umap
#import pickle

import shapely
import scipy
# import fdasrsf

#from matplotlib.lines import Line2D
#from matplotlib.patches import Patch
#import pyarrow.parquet as pq


# import skfda
def cell_features(x):
  x_copy = x.copy()
  mu = x_copy.mean(axis=1)
  for i in range(x_copy.shape[1]):
    x_copy[:,i] = x_copy[:,i] - mu
  scale = np.linalg.norm(x_copy, ord='fro')
  U, sigma, V_T = np.linalg.svd(x_copy)
  major_axis = sigma[0]
  minor_axis = sigma[1]

  eccentricity = np.sqrt( 1 - ( minor_axis ** 2/major_axis ** 2 ))#(major_axis+minor_axis)
 
  area_in_px = shapely.area( shapely.Polygon(x.transpose()) )
  perimeter = shapely.length( shapely.Polygon(x.transpose()) )
  shapely_mp_obj = shapely.MultiPoint( x_copy.T )
  hull_approx_obj = shapely.convex_hull( shapely_mp_obj )
  convexity_score = hull_approx_obj.length/perimeter
  roundness = perimeter ** 2/(4 * np.pi * area_in_px)
          
  return [scale, major_axis, minor_axis, eccentricity, area_in_px,perimeter,roundness,convexity_score]

def preprocess(x):
  mu = x.mean(axis=1)
  for i in range(x.shape[1]):
    x[:,i] = x[:,i] - mu 
  x = x/np.linalg.norm(x, ord='fro')
  return x

def interp_shape_and_preprocess(x_original,num_vertices):
    x = preprocess(x_original)
    x_u = np.unique(x,axis=1)
    x_u, unique_mask = np.unique(x,axis=1, return_index=True)
    points_to_keep = list(set(np.arange(0,x.shape[1])).intersection( set(unique_mask) ))
    x_u = x[:,points_to_keep]
        
    x_coords = np.array(x_u[0,:])
    y_coords = np.array(x_u[1,:])
    tck, u = scipy.interpolate.splprep([x_coords,y_coords], s=0, k=1, per=True)
    xi, yi = scipy.interpolate.splev(np.linspace(0, 1, num_vertices), tck)

    x_interp = np.array([xi,yi])
    x_interp_preprocessed = preprocess(x_interp)

    return x_interp_preprocessed

def reparam_OPA(A, B):
  A_best = A # We start with the guess that A is best aligned to B
  best_error = np.linalg.norm(A - B) # Calculate error between the two matrices
  # Align A to B, consider all possible warping of A
  for shift in range(A.shape[1]):
      gamma = np.roll(np.arange(0, A.shape[1]), shift)
      A_reparam = A[:,gamma] # Try reparameterization of A by composing it with warping gamma
      A_reparam_rot = OPA(A_reparam, B) # Rotationally aligned and reparameterized A with B

      error = np.linalg.norm(A_reparam_rot - B)
      if best_error > error:
          best_error = error
          A_best = A_reparam_rot
  return A_best


def OPA(A, B):
  u, sigma, v_t = np.linalg.svd(B@A.T)
  R = u@v_t
  return R@A

def exp(p, v, theta=None):
  """ Pushes the object on the manifold """
  if not theta:
    theta = np.linalg.norm(v)
  p2 = np.cos(theta)*p + np.sin(theta)*v/np.linalg.norm(v)
  return p2/np.linalg.norm(p2, ord='fro')

def log(p1, p2, reparam=True):
  """ Straightens the pearl beads in a tangent space pointing towards p2 """
  s = p2.shape
  if reparam:
    p2 = reparam_OPA(p2, p1)
  else:
    p2 = OPA(p2, p1)
  p2_v = p2.flatten()#p2.reshape(s[0]*s[1])
  p1_v = p1.flatten()#p1.reshape(s[0]*s[1])
  
  dot = np.clip(np.dot(p1_v, p2_v), -1, 1)
  theta = np.arccos(dot)
  
  # TO calculate the limit and prevent blowup when both vectors are exactly equal
  if theta==0:
    frac = 1
  else:
    frac = (theta/np.sin(theta))

  t = (p2_v - np.dot(p2_v, p1_v)*p1_v)*frac
  return t.reshape(s)

#Assumed that X is a 3-D array of dimensions [# of shapes,2,# of landmarks]
def frechet(X, eta=0.001,use_parallel=False,num_chunks=4,max_iter=100,tol=1e-4,
           store_history=False):
  """ Frechet mean of a set of shapes """
  history = []
  mu = X[0]

  if store_history:
        mu_history = np.zeros((max_iter,X.shape[1],X.shape[2]))
  for i in range(max_iter): # replace that with error margin tolerance algorithm for correctness
    dmu = 0
    
    if use_parallel:
        chunk_size = int(X.shape[0]/num_chunks)
        array_idx_list = np.arange(0,X.shape[0])
        chunks = [array_idx_list[i * chunk_size:(i + 1) * chunk_size] for i in range((len(array_idx_list) + chunk_size - 1) // chunk_size )]

        with Pool(num_chunks) as p:
            dmu_list = p.starmap( log_map_values, zip(repeat(X),repeat(mu),chunks))
            
        for idx in range(len(dmu_list)):
            for x in dmu_list[idx]:
                dmu += x/X.shape[0]
                
    else:
        for j in range(X.shape[0]):
            dmu += log(mu, X[j])/X.shape[0]
           
    prev_mu = mu
    mu = exp(mu, dmu*eta)
    if store_history:
        mu_history[i,] = mu
        
    change_in_mu = theta_value(mu,prev_mu)
    print(change_in_mu,flush=True)
    history.append(change_in_mu)
    if change_in_mu < tol and i > 0:
        break
    
#   plt.plot(history)
    last_itr = i
  if store_history:
    return [mu,mu_history[:last_itr],history]
  else:
    return mu


def PGA(X,mu=None,max_frechet_mean_iter=100,use_parallel=True,frechet_mean_tol=1e-4,
        num_chunks=4,eta=1,store_history=False):
    
    if mu is None:
        if store_history:
            mu_info = frechet(X, eta=eta,max_iter=max_frechet_mean_iter,num_chunks=num_chunks,use_parallel=use_parallel,
                tol=frechet_mean_tol,store_history=True)
            mu = mu_info[0]
        else:
            mu = frechet(X, eta=eta,max_iter=max_frechet_mean_iter,num_chunks=num_chunks,use_parallel=use_parallel,
                tol=frechet_mean_tol,store_history=False)
            
    U = np.zeros(X.shape)
    for i in range(U.shape[0]):
        U[i,:,:] = log(mu, X[i])
        dotprod = np.dot(mu.flatten(), U[i,:,:].flatten())
   
    U_flat = U.reshape(X.shape[0], 2*mu.shape[1])

    P, sigma, R_t = np.linalg.svd(U_flat, full_matrices=False)
    lambdas = (sigma**2)/(X.shape[0]-1)

    if store_history:
        return (P@np.diag(sigma)), lambdas, np.transpose(R_t), mu_info, U_flat
    else:
        return (P@np.diag(sigma)), lambdas, np.transpose(R_t), mu, U_flat

def theta_loop(shapes,ref_idxes,query_idxes):
    num_shapes = len(query_idxes)
    theta_values = np.zeros((len(ref_idxes),len(query_idxes)))
    row_idx = 0
    for ref_idx in ref_idxes:
        col_idx = 0
        for query_idx in query_idxes:
            theta_values[row_idx,col_idx] = theta_value(shapes[ref_idx],shapes[query_idx])
            col_idx += 1
            
        row_idx += 1
        break
        
        
    return theta_values

def theta_value(p1, p2):
    """ Straightens the pearl beads in a tangent space pointing towards p2 """
    s = p2.shape
    p2 = reparam_OPA(p2, p1)
    p2_v = p2.flatten()#p2.reshape(s[0]*s[1])
    p1_v = p1.flatten()#p1.reshape(s[0]*s[1])

    dot = np.clip(np.dot(p1_v, p2_v), -1, 1)
    theta = np.arccos(dot)

    return theta

def log_map_values(shapes,ref_shape,query_idxes):
    num_shapes = len(query_idxes)
    log_map_array = np.zeros((len(query_idxes),shapes.shape[1],shapes.shape[2]))
    ctr = 0
    for idx in query_idxes:
        log_map_array[ctr] = log(ref_shape,shapes[idx])
        ctr += 1
        
    return log_map_array

def theta_loop(shapes,ref_idxes,query_idxes):
    num_shapes = len(query_idxes)
    theta_values = np.zeros((len(ref_idxes),len(query_idxes)))
    row_idx = 0
    for ref_idx in ref_idxes:
        col_idx = 0
        for query_idx in query_idxes:
            theta_values[row_idx,col_idx] = theta_value(shapes[ref_idx],shapes[query_idx])
            col_idx += 1
            
        row_idx += 1
        break
        
        
    return theta_values

        
def plot_triangle(ax, x, edgecolor='black', plot_title=None, plot_type='fill', linewidth=3, point_size=1,
                 label_landmarks=False):
    
    if plot_type == 'fill':
        ax.fill(x[0,:], x[1,:],rasterized=True)
    elif plot_type == 'scatter':
       ax.fill(x[0,:], x[1,:], facecolor='None', edgecolor=edgecolor, linewidth=linewidth,rasterized=True)
       ax.scatter(x[0,:], x[1,:], facecolor='None', edgecolor=edgecolor, linewidth=linewidth, s=point_size,rasterized=True)

    if plot_title is not None:
        ax.set_title(plot_title)
        
    if label_landmarks:
        for i in range(x.shape[1]):
            ax.annotate(str(i), (x[0,i], x[1,i]))

def plot_shape_modes(p,v,mu,lambdas,num_pcs=8,max_sd=3,figsize=(10,5),grid_dims=None,
                     linewidth=1.5,point_size=5,remove_tick_labels=False):

    if grid_dims is None:
        grid_dims = (2,int(num_pcs/2))
        
    fig, axs = plt.subplots(grid_dims[0],grid_dims[1], figsize=figsize,#subplot_kw=dict(box_aspect=1),
                             sharex=True, sharey=True, layout="constrained" )

    # edgecolors = ['green', 'mediumblue', 'darkblue', 'orange', 'red']
    edgecolors = ['gray', 'blue', 'red']#@, 'orange', 'red']

    # lambdas = lambdas * (p.shape[0]-1)
    var_explained = 100 * lambdas/sum(lambdas)

    for j in range(num_pcs):
      scales = np.array([0.000001, +max_sd, -max_sd])
      legend_elements = [Line2D([0], [0], color=color, lw=4, label=str(round(scale, 2))+r'$\sigma$') for color, scale in zip(edgecolors, scales)]

      for i, scale in enumerate(scales):
        if grid_dims[0] > 1:
            if j < np.floor(num_pcs/2):
                ax=axs[0,j]
            else:
                ax=axs[1,j-int(num_pcs/2)]
        else:
            ax=axs[j]

        vec = v[:,j]
        ax.set_title(scale)
        plt.title(np.dot(mu.flatten(), vec))
        theta = scale*np.sqrt(lambdas[j])
        vec = exp(mu, vec.reshape(mu.shape), theta=theta)
        plot_triangle(ax, vec, edgecolor=edgecolors[i],linewidth=linewidth,plot_type='scatter',point_size=point_size)
        #ax.legend(edgecolors, edgecolors)
        ax.legend(handles=legend_elements, loc='upper right')
        ax.set_title('PC{} ({:.2f}%)'.format(j+1,var_explained[j]),pad=0.5)
        # ax.xaxis.label.set_size(15)
        # ax.yaxis.label.set_size(15)
        if remove_tick_labels:
            ax.set_xticklabels([])  # Remove x-axis tick labels
            ax.set_yticklabels([])
            ax.tick_params(axis='x', bottom=False) 
            ax.tick_params(axis='y', left=False) 

        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(linewidth)

        ax.get_legend().remove()

    fig.set_constrained_layout_pads(w_pad=0, h_pad=0.05, hspace=0, wspace=0.1)

def project_onto_dataset(X,mu,v):
    U = np.zeros(X.shape)
    for i in range(U.shape[0]):
        U[i,:,:] = log(mu, X[i])

    U_flat_T = U.reshape(X.shape[0], 2*mu.shape[1])
    return (U_flat_T @ v)

def compute_pre_shape_embedding( boundary_parquet_path,
             pre_shape_output_dir,
             num_vertices_to_sample=30,
         cell_ids_to_analyze=None,
             x_vertex_col="vertex_x", y_vertex_col="vertex_y",
           cell_id_col="cell_id"):
    
    df_shapes = pq.read_table(boundary_parquet_path)

    os.makedirs(pre_shape_output_dir,exist_ok=True)
    
    groups = df_shapes.group_by([cell_id_col]).aggregate(
        [(x_vertex_col,"list"),(y_vertex_col,"list")])

    group_names_array = groups.column(cell_id_col)
    distinct_group_names = group_names_array.to_pylist()
    if cell_ids_to_analyze is None:
        cell_ids_to_analyze = distinct_group_names
        num_groups = groups.num_rows
    else:
        num_groups = len(cell_ids_to_analyze)
    
    num_dimensions = 2
    
    pre_shape_space_embedding_mat = np.zeros((num_groups,num_dimensions,num_vertices_to_sample))
    idx_vals = []
    cell_ids = []
    areas_in_px = np.zeros(num_groups)
    roundness_values = np.zeros(num_groups)
    major_axis_vals = np.zeros(num_groups)
    minor_axis_vals = np.zeros(num_groups)
    perimeters_in_px = np.zeros(num_groups)
    eccentricity_vals = np.zeros(num_groups)
    convexity_vals = np.zeros(num_groups)

    idx = 0
    for group in groups.to_struct_array():
        cell_id = str(group[cell_id_col])
        if cell_id not in cell_ids_to_analyze:
            continue

        if idx % 1000 == 0 and idx > 0:
            print('{}/{}'.format(idx,num_groups),flush=True)
        x_original = np.array( [group['{}_list'.format(x_vertex_col)].values.to_numpy(),
        group['{}_list'.format(y_vertex_col)].values.to_numpy()] )
    
        cell_ids.append( cell_id )
    
        scale, major_axis_vals[idx], minor_axis_vals[idx], eccentricity_vals[idx], areas_in_px[idx],perimeters_in_px[idx],roundness_values[idx], convexity_vals[idx] = cell_features(x_original)
    
        x_interp_preprocessed = interp_shape_and_preprocess(x_original.copy(),num_vertices_to_sample)
        pre_shape_space_embedding_mat[idx] = x_interp_preprocessed
        idx_vals.append( idx )
        idx += 1

            
    idx_df = pd.DataFrame({'cell_id':cell_ids,'numpy_idx':np.array(idx_vals),
                                            'R_idx':np.array(idx_vals)+1,
                       'roundness':roundness_values,
                        'area':areas_in_px,
                         'major_axis':major_axis_vals,
                           'minor_axis':minor_axis_vals,
                      'perimeter':perimeters_in_px,
                      'eccentricity':eccentricity_vals,
		      'convexity':convexity_vals})

    meta_data_file_name = 'Shape_Metadata.csv.gz'
    meta_data_file_path = os.path.join( pre_shape_output_dir, meta_data_file_name )
    idx_df.to_csv( meta_data_file_path, index=False )

    h5_file_name = 'Pre_Shape_Space_Embedding.h5'
    h5file_path = os.path.join( pre_shape_output_dir, h5_file_name )
    
    with h5py.File( h5file_path, 'w'  ) as h5_output:
        h5_output.create_dataset( 'pre_shape_space_embedding', data=pre_shape_space_embedding_mat )

def run_pga( pre_shape_input_dir, pga_output_dir, cell_ids_to_analyze=None, cell_id_col='cell_id',
           max_frechet_mean_iter=1000, eta=1, use_parallel=False,
           num_threads=8,frechet_mean_tol=1e-4):

    os.makedirs(pga_output_dir,exist_ok=True)

    input_h5_file_name = 'Pre_Shape_Space_Embedding.h5'
    pre_shape_file_path = os.path.join( pre_shape_input_dir, input_h5_file_name )
    
    pga_output_file_name = 'PGA_Info.h5' 
    pga_output_file_path = os.path.join( pga_output_dir, pga_output_file_name )

    meta_data_file_name = 'Shape_Metadata.csv.gz'
    meta_data_file_path = os.path.join( pre_shape_input_dir, meta_data_file_name )
    meta_data_df = pd.read_csv( meta_data_file_path )

    if cell_ids_to_analyze is not None:
        mask = meta_data_df[cell_id_col].isin( cell_ids_to_analyze )
        numpy_idxes_to_analyze = meta_data_df.loc[mask,'numpy_idx'].values
    else:
        numpy_idxes_to_analyze = meta_data_df['numpy_idx'].values

    with h5py.File( pre_shape_file_path, 'r' ) as pre_shape_input:
        pre_shape_space_embedding_mat = pre_shape_input['pre_shape_space_embedding'][:]

    with h5py.File( pga_output_file_path, 'w'  ) as h5_output: 
        p, lambdas, v, mu, U_flat = PGA(pre_shape_space_embedding_mat[numpy_idxes_to_analyze,], 
                                        max_frechet_mean_iter=max_frechet_mean_iter,eta=eta,
                                     use_parallel=use_parallel,frechet_mean_tol=frechet_mean_tol,num_chunks=num_threads,store_history=False)
    
        h5_output.create_dataset( 'processed_idxes', data=numpy_idxes_to_analyze)
        h5_output.create_dataset( 'embedding', data=p)
        h5_output.create_dataset( 'variances', data=lambdas)
        h5_output.create_dataset( 'v_matrix', data=v)
        h5_output.create_dataset( 'frechet_mean', data=mu)
        h5_output.create_dataset( 'u_flattened', data=U_flat)
    
def reconstruct_shapes_from_pca(v, mu, lambdas, sds_to_plot, num_pcs=8):

    shapes_df = pd.DataFrame()
    for j in range(num_pcs):
        for i, scale in enumerate(sds_to_plot):
            vec = v[:, j]

            if scale == 0:
                scale = 1e-6
                scale_to_store = 0
            else:
                scale_to_store = scale
                
            # your title logic (kept, but use ax not plt to avoid "current axes" issues)
            theta = scale * np.sqrt(lambdas[j])
            vec2 = exp(mu, vec.reshape(mu.shape), theta=theta)
            df = pd.DataFrame(vec2.T).rename(columns={0:'x',1:'y'})
            df["PC"]     = f"PC{j + 1}"
            df["sd_val"] = float(scale_to_store)
            df["PC"] = df["PC"].to_numpy(dtype=object)
            
            shapes_df = pd.concat( [shapes_df,df], axis=0)

    shapes_df = shapes_df.reset_index(drop=True)
    shapes_df["PC"] = shapes_df["PC"].to_numpy(dtype=object)
    
    return {
        "x":      shapes_df["x"].to_numpy(dtype=float),
        "y":      shapes_df["y"].to_numpy(dtype=float),
        "PC":     shapes_df["PC"].to_numpy(dtype=object),
        "sd_val": shapes_df["sd_val"].to_numpy(dtype=float),
    }

def exp_safe(p, v, theta=None, eps=1e-12):
    """
    Exponential map on the unit sphere/Kendall preshape sphere.

    p: mean shape, shape (2, L)
    v: tangent vector at p, shape (2, L)
    theta: optional geodesic distance. If None, use ||v||.
    """
    v_norm = np.linalg.norm(v)

    if v_norm < eps:
        return p.copy()

    if theta is None:
        theta = v_norm

    p2 = np.cos(theta) * p + np.sin(theta) * v / v_norm
    return p2 / np.linalg.norm(p2, ord='fro')


def reconstruct_shape_from_pca_coords(
    pca_coords,
    v,
    mu,
    lambdas=None,
    coords_are_sds=False,
    return_df=True
):
    """
    Reconstruct a single shape from a vector of PCA coordinates.

    Parameters
    ----------
    pca_coords : array-like, shape (d,)
        PCA coordinates for one shape. If coords_are_sds=False, these are
        ordinary TPCA coordinates/scores. If coords_are_sds=True, these are
        coordinates in standard-deviation units along each PC.

    v : np.ndarray, shape (2L, n_pcs)
        PCA direction matrix returned by PGA. Columns are tangent PC directions.

    mu : np.ndarray, shape (2, L)
        Frechet mean shape.

    lambdas : np.ndarray, optional
        PCA variances. Required if coords_are_sds=True.

    coords_are_sds : bool
        If True, convert SD-unit coordinates to TPCA coordinates by multiplying
        by sqrt(lambda_j).

    return_df : bool
        If True, return a dataframe with x,y coordinates. If False, return
        the reconstructed shape as a (2, L) array.

    Returns
    -------
    pd.DataFrame or np.ndarray
        Reconstructed shape.
    """
    pca_coords = np.asarray(pca_coords).astype(float)
    num_pcs = len(pca_coords)

    if coords_are_sds:
        if lambdas is None:
            raise ValueError("lambdas must be provided when coords_are_sds=True.")
        pca_coords = pca_coords * np.sqrt(lambdas[:num_pcs])

    tangent_vec_flat = v[:, :num_pcs] @ pca_coords
    tangent_vec = tangent_vec_flat.reshape(mu.shape)

    reconstructed_shape = exp_safe(mu, tangent_vec)

    if return_df:
        df = pd.DataFrame(reconstructed_shape.T, columns=["x", "y"])
        df["coord_source"] = "pca_coord_vector"
        return df

    return reconstructed_shape
