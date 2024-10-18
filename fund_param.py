import numpy as np


def x_intensity(lambda_x):
    """
    以下只计算了X射线的连续谱分布,未计算其特征谱分布;特征谱有几个参数未知
    lambda_0 X射线的短波限;
    lambda_x x射线的波长;
    u_t 靶材的质量吸收系数;
    a X射线的出射角;
    z_t 靶材的原子序数;
    t_be X射线光管铍窗厚度;

    """

    def mass_coefficient(E):
        """
        :param E:波长对应的能量值
        :return: 对应波长靶材的质量系数
        """
        l1 = 3.4119
        l2 = 3.1461
        l3 = 3.0038
        l_1 = 3.412
        l_2 = 3.146
        l_3 = 3.004
        ka = 23.22
        c = 16.9573
        k = 23.220
        # 金的2.85存疑
        n = 2.85
        n_1 = 2.73
        n_2 = 2.61439
        # 铜的为瞎编的,原表中是从zn开始不是cu
        n_3 = 2.3554
        m = 0.307
        if E > ka:
            u = c * k * (12.3981 / E) ** n
        elif ka >= np.array(E) > l_1:
            u = c * l1 * (12.3981 / E) ** n_1
        elif l_1 >= E > l_2:
            u = c * l2 * (12.3981 / E) ** n_2
        elif l_2 >= E > l_3:
            u = c * l3 * (12.3981 / E) ** n_3
        elif E > m:
            u = c * 0.66271 * (12.3981 / E) ** 2.6
        else:
            u = None
            print(f"{E}不在计算范围内)")
        return u


    lambda_0 = 1.23981/45
    u_t = mass_coefficient(1.23981/lambda_x)
    a = 45
    z_t = 45
    t_be = 0.0076
    epsilon=((1/lambda_0**1.65)-(1/lambda_x**1.65))*u_t*(1/np.sin(a/180*np.pi))
    C=(1+(1+2.56*10**-3*(z_t**2))**-1)/((1+2.56*10**3*lambda_0*(z_t**-2))*(0.25*epsilon+10**4))
    wab=np.e**(-0.35*lambda_x**2.86*t_be)
    f=(1+C*epsilon)**-2
    X_intensity=2.72*10**-6*z_t*(lambda_x/lambda_0-1)*(1/lambda_x**2)*f*wab
    if lambda_x==(1.2398/20.216):
        lambda_i=1.2398/20.216
        u_0 = lambda_i / lambda_0
        I0 = np.exp(-0.5 * ((u_0 - 1) / (1.17 * u_0 + 3.2)) ** 2) * (
                    3.22 * 10 ** 6 / (9.76 * 10 ** 4 + z_t ** 4) - 0.39) * (u_0 * np.log(u_0) / (u_0 - 1) - 1)
        I_chr = I0 * X_intensity * 50
        return X_intensity+I_chr
    elif lambda_x==(1.2398/22.724):
        lambda_i=1.2398/22.724
        u_0 = lambda_i / lambda_0
        I0 = np.exp(-0.5 * ((u_0 - 1) / (1.17 * u_0 + 3.2)) ** 2) * (
                    3.22 * 10 ** 6 / (9.76 * 10 ** 4 + z_t ** 4) - 0.39) * (u_0 * np.log(u_0) / (u_0 - 1) - 1)
        I_chr = I0 * X_intensity * 50
        return X_intensity+I_chr
    elif lambda_x==(1.2398/2.916):
        lambda_i=1.2398/2.916
        u_0 = lambda_i / lambda_0
        I0 = np.exp(-0.5 * ((u_0 - 1) / (1.17 * u_0 + 3.2)) ** 2) * (
                    3.22 * 10 ** 6 / (9.76 * 10 ** 4 + z_t ** 4) - 0.39) * (u_0 * np.log(u_0) / (u_0 - 1) - 1)
        I_chr = I0 * X_intensity * 50
        return X_intensity+I_chr
    elif lambda_x==(1.2398/2.891):
        lambda_i=1.2398/2.891
        u_0 = lambda_i / lambda_0
        I0 = np.exp(-0.5 * ((u_0 - 1) / (1.17 * u_0 + 3.2)) ** 2) * (
                    3.22 * 10 ** 6 / (9.76 * 10 ** 4 + z_t ** 4) - 0.39) * (u_0 * np.log(u_0) / (u_0 - 1) - 1)
        I_chr = I0 * X_intensity * 50
        return X_intensity+I_chr
    elif lambda_x==(1.2398/2.834):
        lambda_i=1.2398/2.834
        u_0 = lambda_i / lambda_0
        I0 = np.exp(-0.5 * ((u_0 - 1) / (1.17 * u_0 + 3.2)) ** 2) * (
                    3.22 * 10 ** 6 / (9.76 * 10 ** 4 + z_t ** 4) - 0.39) * (u_0 * np.log(u_0) / (u_0 - 1) - 1)
        I_chr = I0 * X_intensity * 50
        return X_intensity+I_chr
    elif lambda_x==(12398/3.144):
        lambda_i=1.2398/3.144
        u_0 = lambda_i / lambda_0
        I0 = np.exp(-0.5 * ((u_0 - 1) / (1.17 * u_0 + 3.2)) ** 2) * (
                    3.22 * 10 ** 6 / (9.76 * 10 ** 4 + z_t ** 4) - 0.39) * (u_0 * np.log(u_0) / (u_0 - 1) - 1)
        I_chr = I0 * X_intensity * 50
        return X_intensity+I_chr
    elif lambda_x==(1.2398/2.697):
        lambda_i=1.2398/2.697
        u_0 = lambda_i / lambda_0
        I0 = np.exp(-0.5 * ((u_0 - 1) / (1.17 * u_0 + 3.2)) ** 2) * (
                    3.22 * 10 ** 6 / (9.76 * 10 ** 4 + z_t ** 4) - 0.39) * (u_0 * np.log(u_0) / (u_0 - 1) - 1)
        I_chr = I0 * X_intensity * 50
        return X_intensity+I_chr
    elif lambda_x==(1.2398/3.002):
        lambda_i=1.2398/3.002
        u_0 = lambda_i / lambda_0
        I0 = np.exp(-0.5 * ((u_0 - 1) / (1.17 * u_0 + 3.2)) ** 2) * (
                    3.22 * 10 ** 6 / (9.76 * 10 ** 4 + z_t ** 4) - 0.39) * (u_0 * np.log(u_0) / (u_0 - 1) - 1)
        I_chr = I0 * X_intensity * 50
        return X_intensity+I_chr
    else:
        return X_intensity
def Spectral_line_score(metal_i):
    """

    :param metal_i: 元素的原子序数
    :return: gka 元素的谱线分数
    """
    gka=[]
    def spectral_line(x):
        return 1/(1+x)
    for item in metal_i:
        if 12<item<20:
            kba=0.001308*item+0.108688
            gka.append(spectral_line(kba))
        elif 20<=item<32:
            kba=0.102+0.0021*item-3*10**-5*item**2+3*10**-7*item**3
            gka.append(spectral_line(kba))
        elif 32<=item<55:
            kba=-0.16838+0.01392*item-1.49634*10**-4*item**2+5.48728*10**-7*item**3
            gka.append(spectral_line(kba))
        # elif item>=74:
        #     kba=0.13868+1.7352*10**-3*item
        #     gka.append(spectral_line(kba))
        elif 92>=item>=55:
            #l系的谱线分数
            fk=0.609-1.619*10**-3*item-0.03248*np.sin(0.161*(item-51)/180*np.pi)
            gka.append(fk)
        else:
            print(f"元素序列为{item}不在谱线分数公式计算范围内")

    return gka
def Absorption_jump_factor(metal_i):
    """

    :param metal_i: 元素的原子序数
    :return: pk 元素的吸收跃迁因子
    """
    # 金的计算结果和一篇博士论文提供的数据不一致
    pk=[]
    for item in metal_i:
        if 11<=item<=50:
            rk=17.54 - 0.6608 * item + 0.01427 * (item ** 2) - 1.1 * 10 ** -4 * (item ** 3)
            jk=(rk-1)/rk
            pk.append(jk)
        elif 51<=item<=83:
            rk=20.03-0.7732*item+0.01159*(item**2)-5.835*10**-5*(item**3)
            jk=(rk-1)/rk
            pk.append(jk)
        else:
            print(f"元素序列为{item}不在吸收跃迁因子公式计算范围内")
    return pk
def mass_absorption_coefficient(E):
    """
    :param E:波长对应的能量值
    :return: 对应波长各元素的质量吸收系数
    c,l1,l2,l3,n_1,n_2,n,kn_3:计算吸收系数的系数
    l_1,l_2,l_3,ka:是各元素的吸收边能量值
    """
    l1 = [14.3528, 3.8058, 1.0961]
    l2 = [13.7336, 3.5237, 0.951]
    l3 = [11.9187, 3.3511, 0.921]
    l_1 =[11.61,3.233,1.022]
    l_2=[11.443,3.15,0.947]
    l_3=[9.743,2.983,0.928]
    ka = [80.725, 25.514, 8.979]
    c = [19.4943, 17.2453, 14.2775]
    k = [80.7249, 25.514, 8.9789]
    # 金的2.85存疑
    n = [2.85, 2.85, 2.85]
    n_1 = [2.65, 2.714, 2.73]
    n_2 = [2.61439, 2.61439, 2.61439]
    # 铜的为瞎编的,原表中是从zn开始不是cu
    n_3 = [2.3554, 2.3554, 2.3554]
    U = []
    for i in range(len(l1)):
        u_value = []
        for item in E:
            if item > ka[i]:
                u = c[i] * k[i] * (12.3981 / item) ** n[i]
            elif ka[i] >= item > l_1[i]:
                u = c[i] * l1[i] * (12.3981 / item) ** n_1[i]
            elif l_1[i] >= item > l_2[i]:
                u = c[i] * l2[i] * (12.3981 / item) ** n_2[i]
            elif l_2[i] >= item > l_3[i]:
                u = c[i] * l3[i] * (12.3981 / item) ** n_3[i]
            elif i==0 and l3[i]>= item>2.206:
                u=c[i]*3.4249*(12.3981/item)**2.575
            else:
                u = None
                print(f"{item}不在计算范围内),{i}")
            u_value.append(u)
        U.append(u_value)
    return U
def sample_mass_coefficient(C_i,lambda_u,a1):
    """

    :param C_i: 各元素的浓度
    :param lambda_u:各元素在临界波长和短波限范围内的质量吸收系数值
    :param a1: 入射角
    :return: 试样s在入射角为a1时不同波长下的有效质量吸收系数
    """
    u_s_lambda=[]
    for i in range(len(lambda_u[0])):
        u_s = (C_i[0]*lambda_u[0][i]/np.sin(a1/180*np.pi)+C_i[1]*lambda_u[1][i]/np.sin(a1/180*np.pi)+
               C_i[2]*lambda_u[2][i]/np.sin(a1/180*np.pi))
        u_s_lambda.append(u_s)
    return u_s_lambda
def sample_lambda_i_mass_coefficient(C_i,lambda_u,a2):
    """

    :param C_i: 各元素的浓度
    :param lambda_u: 各元素在临界波长的质量吸收系数值
    :param a2: 出射角
    :return: 试样s在出射角为a2时不同临界波长下的有效质量吸收系数
    """
    u_s_lambda=[]
    u_s1=0
    u_s2=0
    u_s3=0
    for i in range(len(lambda_u)):
        u_s1 += C_i[i] * lambda_u[i][0] / np.sin(a2/180*np.pi)
        u_s2 += C_i[i] * lambda_u[i][1] / np.sin(a2/180*np.pi)
        u_s3 += C_i[i] * lambda_u[i][2] / np.sin(a2/180*np.pi)
    u_s_lambda.append(u_s1)
    u_s_lambda.append(u_s2)
    u_s_lambda.append(u_s3)
    return u_s_lambda
def obtain_gi(I_intensity,wk,fka,Jk,result):

    k_i = np.array([x * y * z for x, y, z in zip(Jk, wk, fka)])
    deerta_lambda = 0.002
    gi = []
    for i in range(len(lambda_min_max)):
        total_lambda_wi = 0
        for j in range(len(lambda_min_max[i])):
            intensity_i = x_intensity(lambda_min_max[i][j])
            E = 1.23981 / lambda_min_max[i][j]
            u_i = mass_absorption_coefficient([E])
            lambda_wi = u_i[i][0] / result[i][j][i] * intensity_i * deerta_lambda
            total_lambda_wi += lambda_wi
        gi_obtain=I_intensity[i]/k_i[i]/total_lambda_wi
        gi.append(gi_obtain)
    return gi
# c_i 归一化浓度,a1是入射角,a2是出射角,w是立体角,wk是荧光产额,fka是谱线相对强度,jk是吸收跃变因子,lambda_u是各元素在波长λ的吸收系数
# lambda_i_u 是各元素在波长为λ_i的吸收系数，u_s_j是试样s在波长λj处的吸收系数,
def Relative_intensity(c_i,a1,wk,fka,Jk,lambda_u,u,u_s_j,u_s_i):
    """

    :param c_i: 各元素浓度
    :param a1: 入射角
    :param wk: 荧光产额
    :param fka: 谱线分数
    :param Jk: 吸收跃变因子
    :param lambda_u:各元素在不同波长下的质量系数，每个列表代表一个元素在不同波长下的质量系数(au,ag,cu)
    :param u_s_j:不同荧光射线波长下，试样的质量系数,未除过出射角
    :param u_s_i:不同荧光射线波长下,试样的质量系数.已除过出射角
    :return:各元素在其临界激发波长下的强度
    """
    k_i=np.array([x*y*z for x,y,z in zip(Jk,wk,fka)])
    # result_1为au其波长范围内的us,其中len(u)代表不同特征波长,
    # result为一个大列表包含了三个短波限到临界波长范围内的un*,分别为4*3,2*3,7*3,其中4,2,7为短波限和临界波长范围内步进0.2的波长数量
    # beta_ij包含了三中元素其临界波长吸收限内下不同波长时的beta_ij,分别为4*3*3,2*3*3,7*3*3
    lambda_i_u_scaled = [[x / np.sin(a1/180*np.pi) for x in y] for y in u]
    lambda_u_scaled = [[[x / np.sin(a2 / 180 * np.pi) for x in y] for y in z] for z in lambda_u]

    result_1 = np.zeros((len(lambda_u[0][0]),len(u)))
    result_2 = np.zeros((len(lambda_u[1][0]), len(u)))
    result_3 = np.zeros((len(lambda_u[2][0]), len(u)))
    result=[result_1,result_2,result_3]
    beta_1 = np.zeros((len(lambda_u[0][0]),len(u)))
    beta_2 = np.zeros((len(lambda_u[1][0]), len(u)))
    beta_3 = np.zeros((len(lambda_u[2][0]), len(u)))
    beta_ij = [beta_1,beta_2,beta_3]
    for i in range(len(lambda_u_scaled)):
        for j in range(len(lambda_u_scaled[i][0])):
            for w in range(len(lambda_i_u_scaled)):
                col_elements=lambda_u_scaled[i][w][j]
                row_elements = lambda_i_u_scaled[w]
                result[i][j][w] = col_elements + row_elements[i]
    for i in range(len(result)):
        for j in range(len(result[i])):
            for w in range(len(result[i][j])):
                beta_ij[i][j][w] = result[i][j][w] / result[i][j][i] - 1
    gi=obtain_gi(I_intensity=[8286, 8426, 25597],wk=wk,fka=fka,Jk=jk,result=result)
    deerta = np.zeros((len(lambda_u), len(lambda_u)))
    deerta[0, 1] = 1
    deerta[2, 0] = 1
    deerta[2, 1] = 1
    deerta_lambda=0.002
    total_intensity_i=[]
    relative_intensity_i=[]

    lij_1 = np.zeros((len(u), len(lambda_u[0][0])))
    lij_2 = np.zeros((len(u), len(lambda_u[1][0])))
    lij_3 = np.zeros((len(u), len(lambda_u[2][0])))
    l_i_j=[lij_1,lij_2,lij_3]

    for i in range(len(lambda_min_max)):
        for j in range(len(lambda_min_max[i])):
            E = 1.23981 / lambda_min_max[i][j]
            u_i = mass_absorption_coefficient([E])
            u_s = sample_mass_coefficient(c_i, u_i, a1)
            # 要检查lij的最新编码是否正确,以及对deerta的修改
            for k in range(len(u_s_i)):
                if i == k :
                    l_i_j[i][k][j]=0
                else:
                    # 下面的l1[0]是为了避免一个警告，l1计算结果为单个元素的数组,为此加上索引
                    l1 = (1 / np.array(u_s)) * np.log(1 + (np.array(u_s) / u_s_j[k]))
                    l2 = (1 / u_s_i[i]) * np.log(1 + (u_s_i[i] / u_s_j[k]))
                    l_i_j[i][k][j] = l1[0] + l2
                    # print(l1[0],l2)

    for i in range(len(lambda_min_max)):
        intensity = 0
        pure_intensity = 0
        for j in range(len(lambda_min_max[i])):
            deerta_total = 0
            intensity_i = x_intensity(lambda_min_max[i][j])
            E = 1.23981 / lambda_min_max[i][j]
            u_i = mass_absorption_coefficient([E])
            # print(u_i)
            for k in range(len(u_s_i)):
                if deerta[i,k]==0:
                    continue
                else:
                    #  要检查一遍
                    if lambda_min_max[i][j] <= Critical_excitation_wavelength[k]:
                        # print(l_i_j[i][k][j])
                        # print(k_i[k])
                        # print(u_i[k][0])
                        # print(u[k][k])
                        # print(u_i[i][0])
                        # print(l_i_j[i][k][j])
                        # print(c_i[k])
                        deerta_total += 1/2*k_i[k]*u_i[k][0]*u[k][k]/u_i[i][0]*l_i_j[i][k][j]*c_i[k]
                    else:
                        deerta_total += 0
            sum_deerta = 1 + deerta_total
            # print(sum_deerta)
            lambda_wi=u_i[i]/result[i][j][i]*intensity_i*deerta_lambda

            sum_beta = 1 + sum(x*y for x,y in zip(beta_ij[i][j],c_i))
            # if i==0:
            #     print(sum_deerta)
            #     print(sum_beta)
            intensity += lambda_wi * gi[i] * c_i[i] * k_i[i] * sum_deerta / sum_beta
            pure_w_i=u_i[i]/result[i][j][i]*intensity_i*deerta_lambda
            pure_intensity+=pure_w_i*gi[i]*k_i[i]
        # 样本中元素的强度与纯元素的比值
        relative_intensity=intensity/pure_intensity
        print(intensity)
        total_intensity_i.append(intensity)
        relative_intensity_i.append(relative_intensity)
    return relative_intensity_i

def main_fund_param(r_i):
    """

    :param r_i:是样品测试情况下的相对强度
    :return: 满足基本参数法结果
    """

    r_i=[x for x in r_i]
    I_i = [8286, 8426, 25597]
    R_i = [i / I for i, I in zip(r_i, I_i)]
    C_i=[element / sum(R_i) for element in R_i]
    print(R_i)
    # u_s = [sample_mass_coefficient(C_i, x, a1) for x in lambda_u]
    # u特征x射线的质量系数,每个列表代表一个元素在不同荧光X射线波长下的质量系数
    # u_s_i是不同荧光射线波长下,试样的质量系数.已除过出射角
    u_s_i = sample_lambda_i_mass_coefficient(C_i, u, a2)
    # u_s_j是不同荧光射线波长下，试样的质量系数,未除过出射角
    u_s_j = [x * np.sin(a1 / 180 * np.pi) for x in u_s_i]
    relative_intensity_cal = Relative_intensity(c_i=C_i, a1=a1,  wk=wk, fka=fka, Jk=jk,
                                                lambda_u=lambda_u, u=u, u_s_i=u_s_i,
                                                u_s_j=u_s_j)

    C_i_2 = [0 for _ in range(len(r_i))]
    # C_i_1 = [x*y*(1-z)/(x*(y-z)+z*(1-y)) for x,y,z in zip(R_i,C_i,relative_intensity_cal)]
    C_i_1 = [x/y*z for x,y,z in zip(R_i,relative_intensity_cal,C_i)]

    C_i_1 = [item[0] for item in C_i_1]
    print(C_i_1)
    C_i_1 = [item/sum(C_i_1) for item in C_i_1]
    # 没有写第一个迭代时的判断代码，考虑到第一次迭代基本都不可能符合要求，故只在后续迭代中添加了判断语句
    Exman = False
    while not Exman:
        print(C_i_1)
        u_s_i_1 = sample_lambda_i_mass_coefficient(C_i_1, u, a2)
        # u_s_j是不同荧光射线波长下，试样的质量系数,未除过出射角
        u_s_j_1 = [x * np.sin(a1 / 180 * np.pi) for x in u_s_i_1]
        relative_intensity_cal_1 = Relative_intensity(c_i=C_i_1, a1=a1, wk=wk, fka=fka, Jk=jk,
                                                    lambda_u=lambda_u, u=u, u_s_i=u_s_i_1,
                                                    u_s_j=u_s_j_1)
        # print(relative_intensity_cal_1)
        j=0
        #在第一次while循环使用C_i下面的for循环使用C_i和relative_intensity_cal1,其他时候为C_i_1
        for i in range(len(r_i)):
            # print(R_i)
            C_i_2[i]=R_i[i]/relative_intensity_cal_1[i]*C_i_1[i]
            print(C_i_2)
            # print(C_i_1)
            # if abs(C_i_2[i] - C_i_1[i]) < 5 * 10 ** -5:
            #     j += 1
            # if j == (len(r_i)):
            #     Exman = True
        # print(R_i)
        # print(relative_intensity_cal_1)
        # print(C_i_1)
        # print(C_i_2)
        C_i_2 = [element[0] / sum(C_i_2) for element in C_i_2]
        C_i_2 = [x[0] for x in C_i_2]
        if (all(abs(a-b))< 5*10**-5 for a,b in zip(C_i_2,C_i_1)):
            Exman =True
        # print(C_i_2)
        C_i_1 = [x for x in C_i_2]
        # R_i =[x for x in relative_intensity_cal_1]
        # print(R_i)

    return C_i_2
a1=45
a2=45
# 假设我们有一组浓度数据
#样本4
# i_i = [8498,15,501]
# 减去500
# i_i = [8498,15,1]
#样本5
# i_i =[8334,26,613]
# i_i = [8334,26,113]
#样本6
# i_i=[8459,54,668]
# i_i=[8459,54,168]
# 样本7
# i_i=[8420,78,840]
# i_i=[8420,78,340]
# 样本8
# i_i=[8032,209,1079]
# i_i=[8032,209,579]
# 样本9
# i_i=[7699,123,1325]
# i_i=[7699,123,825]
# 样本10
# i_i=[6854,325,1504]
# i_i=[6854,325,1004]
# 样本11
# i_i=[6743,578,1936]
# i_i=[6743,578,1436]
# 样本12
# i_i=[6086,452,2404]
# i_i=[6086,452,1904]
# 样本13
# i_i=[5412,837,3078]
# i_i=[5412,837,2578]
# 样本14
i_i = [5092,503,4899]
# i_i = [5092,503,4399]
# 样本15
# i_i = [4191, 1778, 3943]
# i_i=[4191,1778,3443]
# 样本16
# i_i = [3201, 1415, 6930]
# i_i = [3201, 1415, 6430]
I_i = [8286,8426,25597]
metal_i=[79,47,29]
# 元素i的特征X射线的波长, 该值只与元素的原子序数有关.
E_i=[9.743,22.163,8.046]
lambda_i = [1.23981 / x for x in E_i]
# 不同元素特征波长下au,ag,cu的质量吸收系数
u=mass_absorption_coefficient(E_i)
# 荧光产额，该参数可能存在较大误差,金要用l系的荧光产额k:0.965,l:0.331
wk=[0.331,0.83,0.41]
Critical_excitation_energy=[11.919,25.514,8.979]
Critical_excitation_wavelength=[1.23981/x for x in Critical_excitation_energy]
# v 是光管的激发电压,单位为kv
V=40
lambda_min=1.23981/V
# 各元素的临界激发波长和短波限的范围，步进为0.02
lambda_min_max=[np.arange(lambda_min,x,0.002) for x in Critical_excitation_wavelength]
E_min_max=[1.23981/x for x in lambda_min_max]
# 各元素在临界激发波长和短波限的范围的波长对应的质量吸收系数
# 各元素在不同波长下的质量系数，每个列表代表一个元素在不同波长下的质量系数(au,ag,cu)
lambda_u=[mass_absorption_coefficient(x) for x in E_min_max]
# u_s是在不同元素的临界波长范围内时,其不同波长的试样s的质量系数，已除过入射角
fka=Spectral_line_score(metal_i)
jk=Absorption_jump_factor(metal_i)
Y=main_fund_param(i_i,)
print(Y)


