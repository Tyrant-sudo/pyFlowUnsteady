import time
import os 
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

import pickle
import pandas as pd
import time
from stl import mesh

# 力的角度判断，输入旋翼输出力的大小
# 输入旋翼取截面 0.65 处形状特点，结合对应

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
 
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

# 提取并插值每个截面，计算数值放入 BEMT_program/airfoils 中

import numpy as np
from scipy.interpolate import splev, splprep, interp1d
from scipy.integrate import cumtrapz

def interpolate_airfoil(Q, N, k=3, D=20, resolution=1000):
    '''Q.shape = (,2) Interpolate N points whose concentration is based on curvature. '''
    res, fp, ier, msg = splprep(Q.T, u=None, k=k, s=1e-6, per=0, full_output=1)
    tck, u = res
    # # tck为B样条(权重，坐标值，阶数3），u为参数权重
    uu = np.linspace(u.min(), u.max(), resolution)
    # # 得到插值点坐标
    x, y = splev(uu, tck, der=0)
    dx, dy = splev(uu, tck, der=1)
    ddx, ddy = splev(uu, tck, der=2)
    cv = np.abs(ddx*dy - dx*ddy)/(dx*dx + dy*dy)**1.5 + D

    # 得到1000个点曲率？
    cv_int = cumtrapz(cv, uu, initial=0)
    # # 求曲率+20的积分函数？
    fcv = interp1d(cv_int, uu)
    # # 积分函数连续
    cv_int_samples = np.linspace(0, cv_int.max(), N)
    u_new = fcv(cv_int_samples)
    # 积分函数上的196个点
    x_new, y_new = splev(u_new, tck, der=0)

    xy_new = np.vstack((x_new,y_new)).T

    
    return xy_new


from scipy import spatial

def normal_airfoil(test_prop,save_path, sec = 0,N=200,k=3):
    
    # 格式为 [N,z,x,y]

    section1 = test_prop[sec,:,[1,2]]
    section1 = interpolate_airfoil(section1,N,k,20)

    dist_mat = spatial.distance_matrix(section1.T, section1.T)
    
    # get indices of candidates that are furthest apart
    i, j = np.unravel_index(dist_mat.argmax(), dist_mat.shape)

    posi_xmin = np.min([i,j])
    posi_xmax = np.max([i,j])
    
    xmin      = section1[0][posi_xmin]
    xmax      = section1[0][posi_xmax]
    y_mean    = np.mean(section1[1][posi_xmin])
    section1[1,:] = section1[1,:] - y_mean
    
    y_xmin    = section1[1][posi_xmin]
    y_xmax    = section1[1][posi_xmax] 
    scale     = xmax - xmin
    section1[0,:] = (section1[0,:] - xmin)/scale
    xmin      = section1[0][posi_xmin]
    xmax      = section1[0][posi_xmax]

    angle     = np.arctan((y_xmax - y_xmin)/(xmax - xmin))
    rot_m     = np.array([[np.cos(angle),np.sin(angle)],[- np.sin(angle),np.cos(angle)]])
    section1 = np.dot(rot_m,section1)

    section1[0,:] = 1 - section1[0,:]

    # plt.plot(section1[0],section1[1])
    # plt.scatter([section1[0][i],section1[0][j]],[section1[1][i],section1[1][j]])
    # plt.show()
    
    section1 = section1.T.tolist()
    with open(save_path + str(sec),'w') as f:

        f.write('sec' + str(sec)+'\n')
        for i in section1:
            f.write(str(i[0])+' '+str(i[1])+'\n')

def build_PointsCloud_from_STL(stl_file, output_file, PointsCloud_Direct, Rtip_L, Rhub_L , Num_Zsclices = 14, Num_EachSlices = 200, save_PointsCloudfile = True):
    """
    PointsCloud_Direct means the direct from hub to tip, which will be redirected to [0,0,1]

    Rtip_L, Rhub_L: the position of Rtip and Rhub at the lenth of stl model
    """
    # 读取 STL 文件
    stl_mesh = mesh.Mesh.from_file(stl_file)
    
    # 提取所有顶点
    # 每个三角形由三个顶点组成，每个顶点有三个坐标（x, y, z）
    all_vertices = np.concatenate([stl_mesh.v0, stl_mesh.v1, stl_mesh.v2])

    # 删除重复的顶点
    unique_vertices = np.unique(all_vertices, axis=0)
    # 将点云数据转换为 DataFrame
    
    Rotated_PointsCloud = rotate_point_cloud_to_target(unique_vertices,PointsCloud_Direct, [0, 0, 1])
    L                    = Rotated_PointsCloud[:,2].max() - Rotated_PointsCloud[:,2].min()
    Rtip, Rhub           = Rotated_PointsCloud[:,2].min() + Rtip_L * L, Rotated_PointsCloud[:,2].min() + Rhub_L * L
    z_sclices            = np.linspace(Rhub,Rtip,Num_Zsclices)
    
    Cropped_PointsCloud  = Rotated_PointsCloud[(Rotated_PointsCloud[:,2] >= Rhub) & (Rotated_PointsCloud[:,2] <= Rtip)]
    modydied_PointsCloud = interpolate_and_slice_point_cloud(Cropped_PointsCloud, z_sclices, Num_EachSlices)
    
    if save_PointsCloudfile:
        modydied_PointsCloud.to_csv(output_file, index=False)

    return modydied_PointsCloud

def build_BEM_from_PointsCloud(PointsCloud_file, output_file, root_rR, direct = True, Beta_34 = 20, Num_Blade = 2, Diameter = 0.5):
    
    PointsCloud = pd.read_csv(PointsCloud_file)

    z_min = PointsCloud["z"].min()
    z_max = PointsCloud["z"].max()
    R     = Diameter / 2

    ratio  = (z_max - z_min)/((1 - root_rR) * R)
    root_x = PointsCloud["x"][PointsCloud["z"].argmin()].mean()
    root_y = PointsCloud["y"][PointsCloud["z"].argmin()].mean()
    
    clockwise = -1 if direct is True else 1
    PointsCloud["x"] = (PointsCloud["x"] - root_x) / ratio
    PointsCloud["y"] = (PointsCloud["y"] - root_y) / ratio  * clockwise 
    PointsCloud["z"] = (PointsCloud["z"] - z_min)  / ratio + root_rR * R
    
    Grouped_PC  = PointsCloud.groupby("z")
    
    Title            = "...BEM Propeller..."    
    Num_Sections     = len(Grouped_PC)
    Num_Blade        = Num_Blade
    Diameter         = Diameter * 100
    Beta_34          = Beta_34

    Feather          = 0
    Pre_Cone         = 0
    Center           = [0.0, 0.0 ,0.0]
    Normal           = [-1.0, 0.0, 0.0]

    Radius_R_list    = []
    Chord_R_list     = []
    Twist_list       = []
    # Skew_R_list      = []
    Rake_R_list      = []
    Sweep_list       = []
    t_c_list         = []
    
    Section_list     = []    

    for z_value, group in Grouped_PC:
        group = group.reset_index()
        
        Radius_R_list.append(z_value / R)

        arg1, arg0, chord = find_furthest_points(group)
        Chord_R_list.append(chord / R)
    
        x0,y0      = group["x"][arg0], group["y"][arg0]
        x1,y1      = group["x"][arg1], group["y"][arg1]
        
        # print(x0,y0,x1,y1)
        twist      = np.arctan((y1- y0)/ (x1 - x0 + 1e-10))
        Twist_list.append(twist  * 180/ np.pi)
    
        Rake_R_list.append(y1/R)
        Sweep_list.append(x1/R)
        # Sweep_list:  `sweepdist::Matrix`   : LE sweep distribution (`sweepdist[:, 1] = r/R`, `sweepdist[:, 2] = x/R`
        # Rake_R_list: `heightdist::Matrix`  : LE height distribution (`heightdist[:, 1] = r/R`, `heightdist[:, 2] = z/R`

        SectionPoints,t_c = get_BEMsection(group, arg0, arg1, chord, twist)
        t_c_list.append(t_c)
        Section_list.append(SectionPoints)



    NumSec_34        = int(np.floor(Num_Sections * 3/4))
    twist_drift      = Twist_list[NumSec_34] - Beta_34
    Twist_list       = [_ - twist_drift for _ in Twist_list]

    Skew_R_list      = [0 for _ in range(Num_Sections)]
    Cli_list         = [0 for _ in range(Num_Sections)]
    Axial_list       = [0 for _ in range(Num_Sections)]
    Tangential_list  = [0 for _ in range(Num_Sections)]
    

    data = [[Radius_R_list[_],Chord_R_list[_],Twist_list[_],Rake_R_list[_],Skew_R_list[_],Sweep_list[_],t_c_list[_],Cli_list[_],Axial_list[_],Tangential_list[_]]\
            for _ in range(Num_Sections)]
    
    # 创建文件并写入数据
    with open(output_file, 'w') as file:
        file.write('...BEM Propeller...\n')
        file.write(f'Num_Sections: {Num_Sections}\n')
        file.write(f'Num_Blade: {Num_Blade}\n')
        file.write(f'Diameter: {Diameter:.8f}\n')
        file.write(f'Beta 3/4 (deg): {Beta_34:.8f}\n')
        file.write(f'Feather (deg): {Feather:.8f}\n')
        file.write(f'Pre_Cone (deg): {Pre_Cone:.8f}\n')
        file.write(f'Center: {Center[0]:.8f}, {Center[1]:.8f}, {Center[2]:.8f}\n')
        file.write(f'Normal: {Normal[0]:.8f}, {Normal[1]:.8f}, {Normal[2]:.8f}\n')
        file.write('\nRadius/R, Chord/R, Twist (deg), Rake/R, Skew/R, Sweep, t/c, CLi, Axial, Tangential\n')
        for row in data:
            file.write(', '.join(f'{val:.8f}' for val in row) + '\n')
        file.write('\n')
        for idx, section_data in enumerate(Section_list):
            file.write(f'Section {idx} X, Y\n')
            for point in section_data:
                file.write(', '.join(f'{val:.8f}' for val in point) + '\n')
            file.write('\n')
        
    return PointsCloud

def find_furthest_points(slice_data):
    """
    找到截面上距离最远的两个点。

    :param slice_data: DataFrame或二维numpy数组，包含截面上的点坐标。
    :return: 一个元组，包含距离最远的两个点的坐标。
    """
    from scipy.spatial.distance import pdist, squareform

    if isinstance(slice_data, pd.DataFrame):
        points = slice_data[['x', 'y']].values
    else:
        points = slice_data

    # 计算所有点对之间的距离
    distance_matrix = squareform(pdist(points, 'euclidean'))

    # 找到距离最远的两个点的索引
    furthest_points_idx = np.unravel_index(np.argmax(distance_matrix), distance_matrix.shape)

    return furthest_points_idx[0], furthest_points_idx[1],np.max(distance_matrix)

def get_BEMsection(SectionPoints, arg0, arg1, chord, twist):
    

    SectionPoints["x"] = SectionPoints["x"] - SectionPoints["x"][arg0]
    SectionPoints["y"] = SectionPoints["y"] - SectionPoints["y"][arg0]

    SectionPoints["x"] = SectionPoints["x"] / chord
    SectionPoints["y"] = SectionPoints["y"] / chord
    
    cos_theta = np.cos(-twist)
    sin_theta = np.sin(-twist)

    rotation_matrix = np.array([[cos_theta, -sin_theta], [sin_theta, cos_theta]])

    # 对每个点应用旋转矩阵
    for index, row in SectionPoints.iterrows():
        rotated_point = rotation_matrix.dot(np.array([row['x'], row['y']]))
        SectionPoints.at[index, 'x'] = rotated_point[0]
        SectionPoints.at[index, 'y'] = rotated_point[1]
        
    if SectionPoints["x"].mean() < 0:
        SectionPoints["x"] = - SectionPoints["x"]
        SectionPoints["y"] = - SectionPoints["y"]
    SectionPoints.loc[arg1, 'x'] = 1.0
    SectionPoints.loc[arg1, 'y'] = 0.0

    # 分割为两部分，并按照角度进行排序
    points_below = SectionPoints[SectionPoints['y'] < 0].sort_values(by='x', ascending=False)
    points_above = SectionPoints[SectionPoints['y'] >= 0].sort_values(by='x')
    
    # 合并并重新排序
    reordered_points = pd.concat([points_below, points_above])
    
    Section_array = np.array(reordered_points[["x","y"]])
    t_c           = Section_array[:,1].max() - Section_array[:,1].min()
    Section_array = np.vstack([Section_array[-1],Section_array])
    
    return Section_array,t_c


from scipy.interpolate import griddata

def interpolate_and_slice_point_cloud(point_cloud, z_slices, Num_EachSlices):
    from scipy.spatial import ConvexHull
    # 读取点云数据
    points = point_cloud
    values = point_cloud[:,2]  # 或者其他与点相关的值，如果有的话
    
    # 插值范围（根据实际数据调整）
    x_range = np.linspace(points[:,0].min(), points[:,0].max(), 100)
    y_range = np.linspace(points[:,1].min(), points[:,1].max(), 100)
    x, y = np.meshgrid(x_range, y_range)

    # 创建一个空的DataFrame用于存储所有切片的数据
    all_slices_data = []
    
    # 对每个切片进行插值
    for z in z_slices:
        # 插值
        grid_z = griddata(points, values, (x, y, np.full(x.shape, z)), method='linear')

        # 构建切片数据
        slice_data = pd.DataFrame({'x': x.ravel(), 'y': y.ravel(), 'z': grid_z.ravel()})

        slice_data.dropna(inplace=True)  # 移除插值中产生的 NaN 值
    
        # 提取 x 和 y 值
        if not slice_data.empty:
            slice_points = slice_data[['x', 'y']].values

            # 计算凸包
            hull = ConvexHull(slice_points)
            # 凸包的顶点索引
            hull_indices = hull.vertices
            # 提取凸包的顶点坐标
            Q = slice_points[hull_indices, :]
            
            xy_new = interpolate_airfoil(Q,Num_EachSlices)
            x_new, y_new = xy_new[:,0], xy_new[:,1]
            
            z_new  = np.array([z]).repeat(Num_EachSlices)
            
            slice_data_new = pd.DataFrame({'x': x_new , 'y': y_new , 'z': z_new})
            
        else:
            slice_data_new = slice_data
        # 添加切片数据到总数据帧中
        all_slices_data.append(slice_data_new)

    return  pd.concat(all_slices_data, ignore_index=True)
    
from scipy.spatial.transform import Rotation as R

def rotate_point_cloud_to_target(points, original_direction, target_direction=[0, 0, 1]):
    # 计算旋转轴（叉积）
    axis = np.cross(target_direction,original_direction) 
    # 计算旋转角度
    angle = np.arccos(np.dot(original_direction, target_direction) / 
                      (np.linalg.norm(original_direction) * np.linalg.norm(target_direction)))
    # 创建旋转对象
    rotation = R.from_rotvec(axis * angle / (np.linalg.norm(axis) + 1e-8))
    # 应用旋转到点云
    rotated_points = rotation.apply(points)
    return rotated_points

