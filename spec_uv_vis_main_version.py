from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
from sys import argv
from numpy import power, exp, arange, gradient, array
from pathlib import Path
from time import sleep
class spec_uv_vis:
    def __init__(self): 

        dir=argv[1]
        start=200
        end=400
        fwhm=3226.22
        number_of_points=1000
        self.directory_name=dir

        print('Choose the input file(s):\n')
        print('1 All .out or .log file(s) from the given directory')
        print('2 The .out or .log file(s), from the given directory, selected by myself')
        ans = input()
        if ans=='2':
            print('\nPut any quantity of input file(s) like that: gaus_output_sys1.log gaus_output_sys2.out')
            list_inp_files = [input()]
            for v in list_inp_files:
                self.ls = v.split(' ')
        elif ans=='1':
            r = Path(dir)
            l = []
            tuple_files = (r.glob('*.log'), r.glob('*.out'))
            for c in range(len(tuple_files)):
                for i in tuple_files[c]:
                    i = str(i)
                    i = i.split('/')
                    k = i[len(i)-1]
                    l.append(k)
                self.ls = l

        else:
            print('Answer not supported!!')
            exit()
        
        self.names_str = []
        for names in self.ls:
            names = names[:len(names)-4]
            self.names_str.append(names)
        
        print('\nDefault parameters of Spec uv vis:\n ')
        print(f' 1 start: {start}\n 2 end: {end}\n 3 FHWM: {fwhm}\n 4 number of points: {number_of_points}\n')
        r = input('Do you wanna change any parameter? [Y/N]: ')

        if r in 'yY':
            print('What parameter (s) do you wanna change it? put on like that: 1=300 2=500')
            choose_ls = [input()]

            for values in choose_ls:
                t = values.split(' ')

            ints_ls=[]
            num_ls=[]
            for strings in t:
                ints = int(strings[2:])
                num_ls.append(strings[0])
                ints_ls.append(ints)

            for ind in range(len(num_ls)):
                if num_ls[ind]=='1':
                    start=ints_ls[ind]
                elif num_ls[ind]=='2':
                    end=ints_ls[ind]
                elif num_ls[ind]=='3':
                    fwhm=ints_ls[ind]
                elif num_ls[ind]=='4':
                    number_of_points=ints_ls[ind]
            print('Loading updates...')
            sleep(0.8)
        elif r in 'nN':
            pass
        else:
            print('Answer not supported!!')
            exit()

        self.start = start
        self.end = end
        self.fwhm= fwhm
        self.number_of_points = number_of_points
        try:
            self.get_transitions()
            self.make_spectrum()
            self.make_df()
            self.make_df_normalized()
            self.get_max_wl()   
            self.shell()
            self.contrib_orb()
        except FileNotFoundError:
            print('File(s) not found it!')
            exit()
        print('\nAvailable outputs:\n')
        print(' 1 output data file\n 2 save multiple plots\n 3 show individual plot\n')
        print('What outputs do you want to? put on like that: 1 2 3')
        l_ls = [input()]
        for ll in l_ls:
            u=ll.split(' ') 
        for index in range(len(u)):
            if u[index]=='1':
                print(f'Outputting data file(s) on {self.directory_name}...')
                self.info_file()
                sleep(0.8)
            elif u[index]=='2':
                print('\nDefault parameters of multiple plots:\n\n 1 Experimental csv file: None')
                print(' 2 Lines of experimetal wavelengths of maximal absorption: None')
                print(' 3 Labels of curves: names of input file(s) itself')
                n = input('\nDo you wanna change any these parameters? [Y/N]: ')
                if n in 'nN':
                    self.vline_exp_lambda =None
                    self.labels_to_legend = None
                    self.exp_csv_file  =None
                    self.output_multiple_plots()
                elif n in 'yY':
                    print('What parameters do you wanna change it? put on like that: 1=exp.csv:sep(;) 2=250,340;exp_1,exp2 3=cam-b3lyp,lc-blyp')
                    in_mult = [input()]
                    for inp in in_mult:
                        in_mult_ls  = inp.split(' ')
                    for splited in in_mult_ls:
                        if '1=' in splited:
                            self.exp_csv_file = (splited[2:splited.find(':')], splited[splited.find('(')+1])
                        elif '2=' in splited or '3=' in splited:
                            self.exp_csv_file = None
                    for splited in in_mult_ls:
                        if '2=' in splited:
                            s = splited[splited.find('=')+1:splited.find(';')]
                            lamb_exp_ls_str=s.split(',')
                            s_str = splited[splited.find(';')+1:len(splited)]
                            labels_exp_ls  =s_str.split(',')

                            lamb_exp_ls_float = []
                            for strs in lamb_exp_ls_str:
                                floats = float(strs)
                                lamb_exp_ls_float.append(floats)
                            self.vline_exp_lambda=(lamb_exp_ls_float, labels_exp_ls)

                        elif '1=' in splited or '3=' in splited:
                            self.vline_exp_lambda = None

                    for splited in in_mult_ls:
                        if '3=' in splited:
                            labels_opt  = splited[2:]
                            self.labels_to_legend = labels_opt.split(',')
                        elif '2=' in splited or '1=' in splited:
                            self.labels_to_legend = None

                    self.output_multiple_plots()
                
                print(f'\nOutputting multiple plots file on {self.directory_name}...')
                sleep(0.8)
            elif u[index]=='3':
                print('genereting plot...')
                sleep(0.8)
                self.show_individual_plot()
            else:
                print('Answer not supported!!')
                exit()
        print('Done!')

    def get_transitions(self):

        self.osc_ls_float_final = []
        self.osc_dict_list = [] 
        self.keys_ls=[]
        for file_name in self.ls:
            wl_list = [] 
            osc_float = []
            osc_str_ls = []
            osc_dict = {}
            file = open(f'{self.directory_name}/{file_name}')
            for i in file:
                if ' Excited State' in i:
                    wl_1 = i.split('\n')[0].split(' ')
                    i = wl_1.index('nm')
                    wl = float(wl_1[i-1])
                    wl_list.append(wl)
                    osc_str = wl_1[i+2]
                    osc_str_ls.append(osc_str)
            osc_str_ls_r = list(reversed(osc_str_ls))
            wl_list_r = list(reversed(wl_list)) 
            kprimer = []
            for c1 in range(len(wl_list_r)):
                if wl_list_r[c1] not in kprimer:
                    kprimer.append(wl_list_r[c1])

            kprimer_none = []
            for c in range(len(wl_list_r)):
                if wl_list_r[c] not in kprimer_none:
                    kprimer_none.append(wl_list_r[c])
                else:
                    kprimer_none.append(None)
            for c in range(len(osc_str_ls_r)):
                osc_f = float(osc_str_ls_r[c][2:])
                osc_float.append(osc_f)           
            l1 = []
            for c2 in range(len(osc_float)):
                l2 = []
                for c3 in range(len(osc_float)):
                    if wl_list_r[c2]==wl_list_r[c3]:
                        if kprimer_none[c2]!=None:
                            l2.append(osc_float[c3])
                        else:
                            l2 = None
                l1.append(l2)
            l3 = []
            for a in range(len(l1)):
                if l1[a]!=None:
                    l3.append(l1[a])
            for k in range(len(l3)):
                osc_dict.update({kprimer[k]:l3[k]})
            self.osc_dict_list.append(osc_dict)
            self.osc_ls_float_final.append(osc_float)
            self.keys_ls.append(wl_list_r)
        return self.osc_dict_list
    
    def make_spectrum(self):

        wl_nm=1/self.fwhm
        wl_cm=power(10, 7)/self.fwhm
        A= 2.174*power(10, 8)
        self.final_map_list=[]
        self.epslon_ls = [] 
        for osc_dict in self.osc_dict_list:
            epslon=[]
            final_map={}
            total_map={}
            for wl_ref in osc_dict.keys():
                for fi in osc_dict[wl_ref]:
                    for wl in arange(self.start, self.end, (self.end-self.start)/self.number_of_points): 
                        freq=1/wl
                        freq_ref=1/wl_ref
                        B=fi/wl_cm
                        dif=freq-freq_ref
                        C=power(dif/wl_nm, 2)
                        eps=A*B*exp(-2.7726*C)
                        if wl in total_map.keys(): 
                            x = total_map[wl]
                            x.append(eps)
                            total_map.update({wl:x}) 
                        else:
                            total_map.update({wl:[eps]}) 
            for wl in total_map.keys():
                y=0
                for values in total_map[wl]:
                    y=y+values
                final_map.update({wl:y})
                epslon.append(y)
            self.epslon_ls.append(epslon) 
            self.final_map_list.append(final_map)
        return self.final_map_list 

    def make_df(self):

        self.df_ls=[]
        self.abs_min_ls=[]
        self.abs_max_ls=[]
        for dict in self.final_map_list:
            x = dict.keys()
            y = dict.values()
            d = {}
            d['w'] = x
            d['abs'] = y
            df = DataFrame(d)
            self.df_ls.append(df)
            mini = df['abs'].min()
            self.abs_min_ls.append(mini)
            maxi = df['abs'].max()
            self.abs_max_ls.append(maxi)
            self.abs_max_ls
        return self.df_ls

    def make_df_normalized(self):

        self.df_ls_norm =[]
        for c in range(len(self.df_ls)):
            df_x = self.df_ls[c]['abs']
            self.df_ls[c].insert(column='abs_norm', value=self.df_ls[c]['abs'], loc=2)
            self.df_ls[c]['abs_norm'] = (df_x-df_x.min())/(df_x.max()-df_x.min())
            self.df_ls_norm.append(self.df_ls[c])
        return self.df_ls_norm
    
    def get_max_wl(self):
        
        self.lambda_max_ls = []
        self.y_max_ls = []
        self.y_max_ls_norm = []
        x = arange(self.start, self.end, (self.end-self.start)/self.number_of_points)
        for c in range(len(self.ls)):
            self.x = x
            f_x_ls = self.epslon_ls
            f_x = f_x_ls[c]
            y = array(f_x)
            y_norm = array(self.df_ls_norm[c]['abs_norm'])
            y_max_norm = []
            lambda_max = []
            y_max = []
            first_derivative = (gradient(y, x, edge_order=2), x)
            second_derivative = (gradient(first_derivative[0], first_derivative[1], edge_order=2), first_derivative[1])
            second_derivative = second_derivative[0]
            second_derivative = [round(y, 2) for y in second_derivative]
            for n in range(5, len(second_derivative)-5):
                if second_derivative[n]<-0.5:
                    if second_derivative[n]<second_derivative[n-5]:
                        if second_derivative[n]<second_derivative[n+5]:
                            if second_derivative[n]<second_derivative[n-1]:
                                    if second_derivative[n]<=second_derivative[n+1]:
                                        lambda_max.append(x[n])
                                        y_max.append(y[n])
                                        y_max_norm.append(y_norm[n])
            self.y_max_ls.append(y_max)
            self.lambda_max_ls.append(lambda_max)
            self.y_max_ls_norm.append(y_max_norm)
        return self.lambda_max_ls
    def shell(self):

        shell_ls_list = []
        for c in range(len(self.ls)):
            shell_ls = []
            file = open(f'{self.directory_name}/{self.ls[c]}', 'r')
            r = file.read()
            sp = r.splitlines()
            for n in range(len(sp)):
                if 'RHF' in sp[n]:
                    shell_ls.append(True)
                else:
                    shell_ls.append(False)
            shell_ls_list.append(shell_ls)
        self.shell_final = []
        for count in range(len(self.ls)):
            if True in shell_ls_list[count]:
                self.shell_final.append('closed')
            else:
                self.shell_final.append('open')

    def contrib_orb(self):

        self.contrib_ls = [] 
        self.labels_ls = []
        self.shell_ls= []
        for c in range(len(self.ls)
                       ):
            file_sec = open(f'{self.directory_name}/{self.ls[c]}')
            for lines in file_sec:
                if 'alpha electrons' in lines:
                    line_ls = lines.split(' ')
            i_num_alpha = (line_ls.index('alpha'))-1
            num_alpha = int(line_ls[i_num_alpha])
            lumo_num = num_alpha+1
            ls = []
            if self.shell_final[c]=='closed':
                file = open(f'{self.directory_name}/{self.ls[c]}', 'r')
                r = file.read()
                lines = r.splitlines()
                block = []
                for i in lines:
                    if ' Excited State' in i or '->'in i:
                        block.append(i.strip())
                for j in block:
                    if 'Excited State' in j:
                        l = []
                    else:
                        d = {}
                        k = j.split(' ')
                        if k[2].isnumeric()==True:
                            i_find = int(k.index('->'))
                            key=k[i_find-1]+k[i_find]+k[i_find+1]
                        else:
                            key=k[0]+k[1]
                        d.update({key:f'{power(float(k[len(k)-1]), 2)*200:.2f}%'})
                        l.append(d)
                        if l not in ls:
                            ls.append(l)
                ls_rev = list(reversed(ls))
                self.contrib_ls.append(ls_rev)
                self.labels_ls.append(f'HOMO: {num_alpha}  LUMO: {lumo_num}')   
            else:
                self.s = 'This program just compute contribution of the pairs of orbitals in the electronic excitation in closed Shell systems'
                self.labels_ls.append(self.s)
                self.contrib_ls.append(None)
            
    def info_file(self):

        for c in range(len(self.ls)):
            file= open(f'{self.directory_name}/{self.names_str[c]}_info.txt', 'w')
            file.write(f'  ** Data from file {self.ls[c]} **\n\n')
            a='-'
            file.write(f'{a*55}\n')
            file.write('              Oscillator strengths      wavelength (nm)\n')
            osc_ls=[]
            for keys in self.osc_dict_list[c].keys():
                for values in self.osc_dict_list[c][keys]:
                    osc_ls.append(values)
            if self.shell_final[c]=='closed':
                for i, keys, contribs, osc in zip(range(len(self.keys_ls[c]), 0, -1), self.keys_ls[c], self.contrib_ls[c], osc_ls):
                    if len(str(i))==1:
                        b='0'
                    else:
                        b=''
                    file.write(f'\nExcited State {b}{i}:    {osc:.4f}                  {round(keys, 1)}\n')
                    for c_1 in range(len(contribs)):
                        for keys_lin, values in zip(contribs[c_1].keys(), contribs[c_1].values()):
                            file.write(f'{keys_lin}   {values}\n')
                file.write(f'\n#Labels = {self.labels_ls[c]}')
            else:
                for i, keys, osc in zip(range(len(self.keys_ls[c]), 0, -1), self.keys_ls[c], osc_ls):
                    if len(str(i))==1:
                        b='0'
                    else:
                        b=''
                    file.write(f'\nExcited State {b}{i}:    {osc:.4f}                  {round(keys, 1)}\n')
                file.write(f'\n{self.labels_ls[c]}')
            file.write('\n')
            file.write(f'{a*55}')
            file.write('\n')
            file.write(f'** Data from Gaussian Convolution for FWHM = {self.fwhm} cm⁻¹ **\n\n')
            file.write(f'{a*38}\n')
            file.write('wavelength of maximal absorption (nm)\n')
            for wls in self.lambda_max_ls[c]:
                file.write(f'          {wls:.1f}\n')
            file.write('\n')
            file.write('  w (nm)         Molar abs (L/mol.cm)\n')
            for i in range(len(self.x)):
                file.write(f'  {round(self.x[i], 1)}                 {float(self.epslon_ls[c][i]):.4f}\n')   
            file.write(f'{a*38}')   
    def show_individual_plot(self):

        x_vline_ls_final = []
        y_vline_ls_final =[]
        ls_final_osc = []
        ls_final_wl = []
        label_ls_final = []
        for c in range(len(self.df_ls)):
            plt.figure(figsize=[6,6], num=self.names_str[c])
            x = self.df_ls[c]['w']
            osc_max = max(self.osc_ls_float_final[c])
            x_vline_ls = []
            y_vline_ls = []
            for keys in self.keys_ls[c]:
                x_vline_ls.append(keys)
            x_vline_ls_final.append(x_vline_ls)
            y_vline_ls_final.append(y_vline_ls)
            for osc in self.osc_ls_float_final[c]:
                osc_r = osc/osc_max
                y_vline = self.abs_max_ls[c]*osc_r
                y_vline=osc_r
                y_vline_ls.append(y_vline)
            y = self.df_ls_norm[c]['abs_norm']
            add = {'family':'serif','color':'black','size':15}
            plt.ylabel('Relative intensity', fontdict=add)
            plt.ylim(0.0, 1.1)
            for count_n in range(len(self.y_max_ls_norm[c])):
                par = self.y_max_ls_norm[c][count_n]
                plt.annotate(xy=(self.lambda_max_ls[c][count_n]+4, par+0.073), 
                                text=f'{self.lambda_max_ls[c][count_n]:.1f} nm', color='blue', 
                                xytext=(self.lambda_max_ls[c][count_n]+4, par+0.073),
                                 bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.2), size=8.5)
            for j in range(len(self.osc_dict_list[c])):
                plt.vlines(x=x_vline_ls_final[c][j], ymin=0.0, ymax=y_vline_ls_final[c][j], colors='black')
            ls_osc = []
            ls_wl = []
            label_ls = []
            for p in range(len(y_vline_ls_final[c])):
                if y_vline_ls_final[c][p]>0.03:
                    ls_osc.append(y_vline_ls_final[c][p])
                    ls_wl.append(x_vline_ls_final[c][p])
                    label_ls.append(len(y_vline_ls_final[c])-p)
            ls_final_osc.append(ls_osc)
            ls_final_wl.append(ls_wl)
            label_ls_final.append(label_ls)
            for k in range(len(ls_final_osc[c])):
                plt.annotate(xy=(ls_final_wl[c][k]-1.45, ls_final_osc[c][k]+0.01), text=f'{label_ls_final[c][k]}', color='red')
            plt.plot(x, y, color='black')
            plt.xlabel('wavelength (nm)', fontdict=add)
            plt.grid(alpha=0.4, linestyle='--')
            plt.xlim(self.start, self.end) 
            plt.show()
            
    def output_multiple_plots(self):               

        colors_new = []
        colors = ["#020E07",'#EC2504', "#12D10C", 
                    "#09C4F2", "#C507CC", "#1F08EA", 
                    "#FAE606", "#F98501", '#CF95D7']
        for c in range(len(self.ls)):
            colors_new.append(colors[c])
            curve_colors=colors_new
        if self.labels_to_legend==None:
            self.labels_to_legend=self.names_str
        plt.figure(figsize=[7,7], num='plots_spec_uv_vis')
        ax = plt.subplot()
        for c in range(len(self.df_ls)):
            x = self.df_ls[c]['w']
            y = self.df_ls_norm[c]['abs_norm']
            ax.plot(x, y, color=curve_colors[c], 
                    label=self.labels_to_legend[c]
                    )  
        if self.exp_csv_file!=None:
            name = self.exp_csv_file[0][:len(self.exp_csv_file[0])-4]
            df = read_csv(f'{self.directory_name}/{self.exp_csv_file[0]}', sep=self.exp_csv_file[1])
            x = df[df.columns[0]]
            y = df[df.columns[1]]
            y = (y-y.min())/(y.max()-y.min())
            ax.plot(x, y, linestyle='--', label=name, color='purple')
        add={'family':'serif','color':'black','size':14}
        plt.ylabel(ylabel='Relative intensity', fontdict=add)
        plt.ylim(0.0, 1.005)
        ax.grid(alpha=0.4, linestyle='--')
        plt.xlabel(xlabel='wavelength (nm)', fontdict=add)
        plt.xlim(self.start, self.end)
        if self.vline_exp_lambda!=None:
            for i in range(len(self.vline_exp_lambda[0])):
                plt.vlines(x=self.vline_exp_lambda[0][i], ymin=0.0, ymax=1.005, linestyles='--', colors='blue', label=self.vline_exp_lambda[1][i])
        if len(self.ls)>=3:
            ax.legend(loc="lower left", bbox_to_anchor=(0.03, 1.00, 0.90, 0.15), mode='expand', ncol=3, frameon=False)
        else:
            ax.legend(loc='upper center', mode='expand', bbox_to_anchor=(0.16, 0.97, 0.65, 0.15), ncol=2, frameon=False)
        plt.savefig(fname=f'{self.directory_name}/plots_spec_uv_vis.png', format='png')
        plt.close()

if __name__== '__main__':
    spec_uv_vis()

