import time
import os 
import numpy as np
import csv
import ast
import sys
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import pickle
import pandas as pd
import time
import math

# import open3d as o3d
from scipy.spatial import ConvexHull
import imageio.v2 as imageio
from scipy.interpolate import splev, splprep, interp1d
from scipy.integrate import cumtrapz

data_root_path   = 'F:\\graduate_student\\T2_GANpropeller\\test2\\1_model\\grid_mesh\\'
def get_group(name_list,save_path,pic_width,pic_heigh,figsize=(20,6)):
    column,row = len(name_list),len(name_list[0])
    
    fig = plt.figure(figsize=figsize)

    for i in range(column):
        for j in range(row):

            posi = j + i * row + 1
          
            ax = fig.add_subplot(column,row,posi)
            ax.axis('off')
            im   = img.imread(name_list[i][j])
            s = im.shape
            left  = (s[0] - pic_width)//2
            right = (s[0] + pic_width)//2
            below = (s[1] - pic_heigh)//2 
            upper = (s[1] + pic_heigh)//2

            im   = im[left:right ,below:upper,:]
            ax.imshow(im)
    plt.savefig(save_path)
    plt.show()
    plt.close()
    
def get_gif(path,save_name):
    all_files = os.listdir(path)

    # 仅保留以".png"结尾的文件
    png_files = [file for file in all_files if file.endswith('.png')]
    images = []

    for i in range(len(png_files)):
        images.append(imageio.imread(path + png_files[i]))

    imageio.mimsave(path + save_name, images, duration=0.1)

def data_root():
    """
    name1 = 'mesh_data_modified.npy' (44467,18,68,3)\n
    name2 = 'mesh_data_test.npy' (4452,18,68,3)
    """

    return data_root_path

def fluent(data):
    new_x = np.linspace(0,1,200)
    x     = np.linspace(0,1,data.shape[0])
    f = interpolate.interp1d(x,data,'cubic')
    new_data = f(new_x)
    return new_data
def interpolate_sec(Q, N=200, k=3, D=20, resolution=1000):
    ''' Interpolate N points whose concentration is based on curvature. '''
    Q = np.concatenate((Q,Q[0][None,:]),0)
    res, fp, ier, msg = splprep(Q.T, u=None, k=k, s=1e-6, per=0, full_output=1)
    tck, u = res
    uu = np.linspace(u.min(), u.max(), N)
    x_new, y_new = splev(uu, tck, der=0)
    xy_new = np.concatenate((x_new[:,None], y_new[:,None]),1)
    return xy_new

def create_line(x,y,z,ax,chord,R):
    lenth = x.shape[0]
   
    z_mean = np.mean(z)
    y_mean = np.mean(y)
    x_mean = np.mean(x)

    center_z = z_mean - R*3/3
    center_y = y_mean
    center_x = x_mean
    center = np.array([[center_x,center_y,center_z]]).repeat(lenth,0)
    
    theta_left = np.arcsin(chord/R)
    theta_left = 0.002
    theta = np.linspace(-theta_left,theta_left,lenth)

    dis_x = np.zeros(x.shape[0])
    displace = np.array([dis_x,R*np.sin(theta),R*np.cos(theta)]).T

    line  = center +  displace
    
    line = line[::3,:]
    # ax.scatter3D(line[:,0],line[:,1],line[:,2],s = 1,c = 'cyan')

    return line


def draw_geom_line(foil_data,name,pic_path,elev = -65,azim = 10):

    pic_name = pic_path +'/'+ str(name) + '.png'
    fig = plt.figure(figsize=(45, 30))
    ax = plt.axes(projection="3d")
    ax.set_xlim(0,1)
    ax.set_ylim(-0.5,0.5)
    ax.set_zlim(-0.5,0.5)

    ax.view_init(elev=elev, azim=azim)

    outline = []

    lines = []
    for i in range(foil_data.shape[0]):

        x = foil_data[i,:,0]
        y = foil_data[i,:,1]
        z = foil_data[i,:,2]

        y_left  = np.argmin(y)
        y_right = np.argmax(y)

        tmp     = [[x[y_left],y[y_left],z[y_left]],[x[y_right],y[y_right],z[y_right]]]
        p_range = outline.append(tmp)

        ax.scatter3D(x,y,z,s = 5,c = 'tomato',alpha = 0.5)
        line = create_line(x,y,z,ax,np.abs(y[y_right] - y[y_left]),100)
        lines.append(line)
        ax.plot3D(x,y,z, c = 'gray',alpha = 0.5)
        
    lines = np.array(lines)
    outline = np.array(outline)
    x_left = fluent(outline[:,0,0])
    y_left = fluent(outline[:,0,1])
    z_left = fluent(outline[:,0,2])
    
    x_right = fluent(outline[:,1,0])
    y_right = fluent(outline[:,1,1])
    z_right = fluent(outline[:,1,2])

    ax.plot3D(x_left,y_left,z_left, c = 'k',alpha = 0.5)
    ax.plot3D(x_right,y_right,z_right, c = 'k',alpha = 0.5)
    
    plt.axis('off')
    plt.savefig(pic_name)
    # plt.show()
    plt.close()
    return lines

def easy_draw(foil_data,save_name,elevation_angle=0, azimuth_angle = 90 ):
    # p = foil_data
    x = foil_data[:,0]
    y = foil_data[:,1]
    z = foil_data[:,2]

    fig = plt.figure(figsize=(10, 7))
    ax = plt.axes(projection="3d")
    max_range = np.array([x.max()-x.min(), 
                          y.max()-y.min(), 
                          z.max()-z.min()]).max() / 2.0

    # 找到中心点
    mid_x = (x.max()+x.min()) * 0.5
    mid_y = (y.max()+y.min()) * 0.5
    mid_z = (z.max()+z.min()) * 0.5

    # 设置统一的轴比例
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.zaxis.set_visible(False)

    ax.scatter3D(x,y,z,s=0.5)
    ax.view_init(elev=elevation_angle, azim=azimuth_angle)

    # plt.axis('off')
    plt.show()
    plt.savefig(save_name)
    plt.close()
    return 


def draw_surface(foil_data,name,pic_path):
    # p = foil_data
    x = foil_data[:,:,0]
    y = foil_data[:,:,1]
    z = foil_data[:,:,2]
    pic_name = pic_path +'/'+ str(name) + '.png'

    fig = plt.figure(figsize=(10, 7))
    ax = plt.axes(projection="3d")
    ax.set_xlim(0,1)
    ax.set_ylim(-0.5,0.5)
    ax.set_zlim(-0.5,0.5)

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.zaxis.set_visible(False)

    ax.plot_surface(x,y,z,cmap= 'ocean')

    plt.axis('off')
    plt.savefig(pic_name)
    plt.close()


def distance(p1, p2):
    """
    计算两点之间的距离
    """
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

def angle(p1, p2):
    """
    计算极角
    """
    return math.atan2(p2[1]-p1[1], p2[0]-p1[0])

def clockwise_sort(points):
    """
    按照从近到远且顺时针的顺序排序点集
    """
    # 计算中心点
    center_x = sum([p[0] for p in points]) / len(points)
    center_y = sum([p[1] for p in points]) / len(points)
    center = (center_x, center_y)

    # 计算每个点到中心点的距离和极角
    polar_points = [(p, distance(p, center), angle(center, p)) for p in points]

    # 按照极角排序
    polar_points = sorted(polar_points, key=lambda x: x[2])

    # 如果存在多个点共线的情况，按照距离排序
    for i in range(len(polar_points)-1):
        p1, dist1, angle1 = polar_points[i]
        p2, dist2, angle2 = polar_points[i+1]
        if angle1 == angle2 and dist1 > dist2:
            polar_points[i], polar_points[i+1] = polar_points[i+1], polar_points[i]

    # 返回排序后的点集
    a = [list(p[0]) for p in polar_points]
    
    return np.array(a)

"""
def get_value_center(data):



    pcd = o3d.geometry.PointCloud()
    
    pcd.points = o3d.utility.Vector3dVector(data)


    tetra_mesh, pt_map = o3d.geometry.TetraMesh.create_from_point_cloud(pcd)
    alpha = 0.05

    mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(
        pcd, alpha, tetra_mesh, pt_map)
    mesh.compute_vertex_normals()
    is_watertight = mesh.is_watertight()
    is_inter      = mesh.is_self_intersecting()

    try:
        volume = mesh.get_volume()
        center = mesh.get_center()
    except:
        return is_inter,False,-1,-1
    return is_inter,is_watertight,volume,center

"""""
"""
def get_section(data,posi = [0.5]):


    pcd = o3d.geometry.PointCloud()
    
    pcd.points = o3d.utility.Vector3dVector(data)


    tetra_mesh, pt_map = o3d.geometry.TetraMesh.create_from_point_cloud(pcd)
    alpha = 0.05

    mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(
        pcd, alpha, tetra_mesh, pt_map)
    mesh.compute_vertex_normals()
    mesh0 = mesh.sample_points_uniformly(500000)
    

    points= np.asarray(mesh0.points)
    
    for i,_ in enumerate(posi):
        po = np.argwhere(np.abs(points[:,0]-_) < 0.001)
        if po.shape[0] <2:
            continue
        p = points[po].reshape(-1,3)
        
        x = fluent(p[:,0])[:,None]
        x = np.ones_like(x) * _

        yz = clockwise_sort(p[:,[1,2]])

        yz = interpolate_sec(yz)
        
        p1 = np.concatenate((x,yz),1)
            
        try:
            result = np.concatenate((result,p1[None,:]),0)
        except:
            result = p1[None,:]

    return result

"""    
    
from scipy import spatial
def draw_section(foil_data,sec_po = [0.3,0.5],name = '',pic_path = '',elev=-20, azim=15,scale = False,displace = 0.5,rota = True):
    
    def farthest_points(points):
        
        pts = points
        
        # two points which are fruthest apart will occur as vertices of the convex hull
        candidates = pts[spatial.ConvexHull(pts).vertices]
        
        # get distances between each pair of candidate points
        dist_mat = spatial.distance_matrix(candidates,candidates)
        
        # get indices of candidates that are furthest apart
        i,j = np.unravel_index(dist_mat.argmax(),dist_mat.shape)
        
        a,b = np.argwhere(points == candidates[i])[0][0],np.argwhere(points==candidates[j])[0][0]
        
        if points[a][0] < points[b][0]:
            return a,b
        else:
            return b,a

    def norm_sec(section,scale,displace):
        
        x  = section[:,0][:,None]
        y  = section[:,1]
        z  = section[:,2]
        points = section[:,[1,2]]
        
        p0,p1 = farthest_points(points)
        # p0 = np.argmin(y)
        # p1 = np.argmax(y)
        points[:,1] = points[:,1] - points[p0,1] 
        points[:,0] = points[:,0] - points[p0,0]
        if scale == True:
            scale = 1/np.sqrt((z[p1]-z[p0])**2 + (y[p1]-y[p0])**2)
        else:
            scale = 1.5
        
        theta = -math.atan2(z[p1]-z[p0],y[p1]-y[p0])
        rot_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                       [np.sin(theta), np.cos(theta)]])
        
        scale_matrix = np.array([[scale, 0], [0, scale]])

        if rota:
            rotated_points = np.dot(rot_matrix, points.T).T
        else:
            rotated_points = points
        scaled_points = np.dot(scale_matrix, rotated_points.T).T

        points = np.concatenate((x,scaled_points),1)
        points[:,1] = points[:,1] - displace
        
        return points
    
    rlt = []
    foil_data1 = get_section(foil_data,sec_po)
    for _ in foil_data1:
        
        standard = norm_sec(_,scale,displace)
        rlt.append(standard)
                
    
    rlt = np.array(rlt)

    pic_name = pic_path +'/'+ str(name) + '.png'
    
    foil_data = foil_data.reshape(-1,68,3)
    x  = foil_data[:,:,0]
    y  = foil_data[:,:,1]
    z  = foil_data[:,:,2]

    
    fig = plt.figure(figsize=(10, 7))
    ax = plt.axes(projection="3d")
    ax.set_xlim(0,1)
    ax.set_ylim(-0.5,0.5)
    ax.set_zlim(-0.5,0.5)
    ax.view_init(elev=elev, azim=azim)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.zaxis.set_visible(False)

    for _ in rlt:
        x1 = _[:,0]
        y1 = _[:,1]
        z1 = _[:,2]
        ax.plot3D(x1,y1,z1,c= 'navy')
    ax.plot_surface(x,y,z,cmap= 'gist_yarg',alpha = 0.7)
    plt.axis('off')
    plt.savefig(pic_name)
    # plt.show()
    plt.close()
    return rlt

def draw_list(crlist,ax,color = 'b'):
    
    hull = spatial.ConvexHull(crlist, qhull_options="QJ")
    poly = plt.Polygon(crlist[hull.vertices, :],alpha = 0.1, color = color)
    ax.add_patch(poly)

def get_data(g):
    
    crlist,betalist = g[:,[0,1]],g[:,[2,3]]
    
    return crlist,betalist

import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D, art3d

from stl import mesh

def plot_stl(filename,save_path,azimuth_angle=0, elevation_angle=0):
    # 读取 STL 文件
    your_mesh = mesh.Mesh.from_file(filename)

    # 创建一个新的绘图
    figure = plt.figure(figsize=(10,10))
    axes = figure.add_subplot(111, projection='3d')

    # 添加 STL 文件的顶点到绘图中
    axes.add_collection3d(art3d.Poly3DCollection(your_mesh.vectors))

    # 计算每个轴的范围
    max_range = np.array([your_mesh.x.max()-your_mesh.x.min(), 
                          your_mesh.y.max()-your_mesh.y.min(), 
                          your_mesh.z.max()-your_mesh.z.min()]).max() / 2.0

    # 找到中心点
    mid_x = (your_mesh.x.max()+your_mesh.x.min()) * 0.5
    mid_y = (your_mesh.y.max()+your_mesh.y.min()) * 0.5
    mid_z = (your_mesh.z.max()+your_mesh.z.min()) * 0.5

    # 设置统一的轴比例
    axes.set_xlim(mid_x - max_range, mid_x + max_range)
    axes.set_ylim(mid_y - max_range, mid_y + max_range)
    axes.set_zlim(mid_z - max_range, mid_z + max_range)
    axes.view_init(elev=elevation_angle, azim=azimuth_angle)
    # 显示绘图
    figure.savefig(save_path)

def draw_BEMgeom():
    
    1


if __name__ == "__main__":

    plt_path = '../0_model/stl/'
    plt_name = 'Exp_16inch.stl'

    file_name = plt_path + plt_name
    save_path = '../2_pic/test/'
    
    a_list = [_ for _ in range(0,180,5)]
    e_list = [_ for _ in range(0,180,5)]
    
    # [0,45] and [90,90]

    for a in a_list:
        for e in e_list:
            save_name = save_path + f'{a}_{e}.png'
            plot_stl(file_name,save_name,a,e)