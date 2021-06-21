#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void double_me(int* x) {
  // Doubles the value at the memory location pointed to by x
  *x = *x + *x;
}

void for_test(int* iter_num, int* vec){
	int i;
	for(i =0; i < iter_num[0]; i++){
		vec[i] = i;
	}
}



void vec_subtract(double* vec, double* a, int* n){
	int i;
	for(i =0; i < n[0]; i++){
		vec[i] = vec[i] - a[i];
	}	
}

void vec_scalar_subtract(double* vec, double* a, int* n){
	int i;
	for(i =0; i < n[0]; i++){
		vec[i] = vec[i] - a[0];
	}	
}

void vec_abs(double* r, int* n){
	int i; 
	for(i = 0; i < n[0]; i++){
		r[i] = fabs(r[i]);
	}
}


void rq_loss(double* r, double* tau, int* n){
	int i; 
	for(i = 0; i < n[0]; i++){
		r[i] = (fabs(r[i])+(2 * *tau - 1) * r[i])/2;
	}
}


void huber_loss(double* r, double* gamma, int* n){
	int i;
	double r_abs;
	for(i=0; i < n[0]; i++){
		r_abs = fabs(r[i]);
		if(r_abs <= *gamma){
			r[i] = r[i]*r[i]/(2 * *gamma);
		} else{
			r[i] = r_abs - *gamma/2; 
		}
	}
}



void rq_huber(double* r, double* tau, double* gamma, int* n){
	int i;
	double r_abs; 
	for(i = 0; i < n[0]; i++){
		r_abs = fabs(r[i]);
		if(r_abs <= *gamma){
			r[i] = r[i]*r[i]/(2 * *gamma) + (2 * *tau - 1)*r[i]/2;
		} else{
			r[i] = r_abs - *gamma/2 + (2 * *tau - 1)*r[i]/2; 
		}
	}
}


void vec_multiply(double* vec, double* a, int* n){
	int i;
	for(i =0; i < n[0]; i++){
		vec[i] = vec[i] * a[i];
	}
}

void rq_huber_deriv(double* r, double* tau, double* gamma, int* n){	
	int i; 
	double r_abs; 
	for(i = 0; i < n[0]; i++){
		r_abs = fabs(r[i]);
		if(r_abs <= *gamma){
			r[i] = r[i]/ *gamma + (2 * *tau -1)/2;
		} else{
			r[i] = copysign(1.0, r[i]) + (2 * *tau - 1)/2;
		}
	}	
}


void rq_tanh_deriv(double* r, double* tau, double* gamma, int* n){
	int i; 
	for(i=0; i < n[0]; i++){
		r[i] = (tanh(r[i]/ *gamma) + 2 * *tau - 1)/2;
	}
}

void m_v_mult_avg(double* m, double* v, int* n, int* p, double* ret){
	int i,j;
	for(j=0; j < p[0]; j++){
		for(i = 0; i < n[0]; i++){
			ret[j] = ret[j] + m[j*n[0]+i]*v[i]/n[0];	
		}
	}
}


void neg_gradient(double* r, double* weights, double* tau, double* gamma, double* x, int* n, int p, int* apprx, int x_start_pos, double* ret){
	int i, j; 
	double grad_val;
	//printf("n and p %d %d \n", n[0], p);
	for(j = 0; j < p; j++){
		ret[j] = 0; 
		for(i = 0; i < n[0]; i++){
			grad_val = 0; 
			//printf("Working on j i %d %d \n", i, j);
			//ret[j] = ret[j]+i;
			if(apprx[0] == 1){
				//printf("using huber");
				if(fabs(r[i]) <= gamma[0]){
					grad_val = (r[i]/ *gamma + (2 * *tau -1))/2;
				} else{
					grad_val = (copysign(1.0, r[i]) + (2 * *tau - 1))/2;
				}
			} else{
				//printf("using tanh");
				grad_val = (tanh(r[i]/ *gamma) + 2 * *tau - 1)/2;
			}
			ret[j] = ret[j]+grad_val * weights[i] * x[x_start_pos+j*n[0]+i]/n[0];	
			//printf("Contribution at predictor sample value %d %d %f \n", j, i, grad_val * weights[i] * x[x_start_pos+j*n[0]+i]/n[0]);
			//printf("Group sample grad_val %d %d %f \n", j, i, grad_val);
			//printf("x pos value is %d %f \n", j*n[0]+i, x[x_start_pos+j*n[0]+i]);
		}
		//printf("In neg_gradient function pos value %d %f \n", j, ret[j]);
	}
}

void neg_gradient_r_test(double* r, double* weights, double* tau, double* gamma, double* x, int* n, int* p, int* apprx, int* x_start_pos, double* ret){
	neg_gradient(r,weights, tau, gamma, x, n, p[0],apprx,x_start_pos[0],ret);
}


//void print_array(int* vec, int* n){
	//int i; 
	//for(i = 0; i < *n; i++){
		//printf("pos val %d %d \n", i, vec[i]);
	//}
//}

//void print_dbl_array(double* vec, int n){
	//int i; 
	//for(i = 0; i < n; i++){
		//printf("pos val %d %f \n", i, vec[i]);
	//}
//}

void l2norm(double *r, int n, double* norm){
	int i;
	//printf("in l2 norm printing array \n");
	//print_dbl_array(r,n);
	for(i = 0; i < n; i ++){
		//printf("in l2norm pos value_added %d %f \n", i, r[i]*r[i]);
		norm[0] = norm[0] + r[i] * r[i];
		//printf("new norm value is %f \n", norm[0]);
	}
	norm[0] = sqrt(norm[0]);
}



int inarray(int* val, int *arr, int* size){
    int i;
    for (i=0; i < size[0]; i++) {
        if (arr[i] == *val){
            return 1;
		}
    }
    return 0;
}

// calculates group norms correctly, but we don't need to be calculating all of them. I think we should subset it at the R level. 
// Need to verify that this gets done correctly. Maybe..., who cares. 

// void kkt_check(double* r, double* weights, double* inactive_w_lambda, double* gamma, double* tau, int* inactive_g,
				// int* inactive_ind, double* lambda, double* inactive_x, int n, int inactive_p, int* inactive_g_num, int* apprx, int* fail_idx, int* kkt_fail){
	// double grad_ret[inactive_p[0]];
	// double group_norms[inactive_g_num[0]];
	// neg_gradient(r, weights, tau, gamma, inactive_x, n, inactive_p, apprx, grad_ret);
	// /*int i;
	// for(i=0; i < inactive_p[0]; i++){
		// printf("%f\n",grad_ret[i]);
	// }
	// //for(i=0; i < )
	// /*int idx;
	// printf("Before group norms\n");	
	// printf("%d\n", p[0]);*/
	// int i, j;
	// double compare_val;
	// for(i=0; i < inactive_g_num[0]; i++){
		
		// for(j = 0; j < inactive_p[0]; j++){
			// if(inactive_ind[j]==inactive_g[i]){
				// /*printf("predictor: ");
				// printf("%d\n",j);
				// printf(" is part of group ");
				// printf("%d\n", inactive_g[i]);
				// printf("with an added value of ");
				// printf("%f\n", grad_ret[j]*grad_ret[j]);*/
				// group_norms[i] = group_norms[i] + grad_ret[j]*grad_ret[j];
			// }
		// }
		// /*printf("for group");
		// printf("%d\n", i);
		// printf("group norm value of ");
		// printf("%f\n", group_norms[i]);
		// printf("with a lambda weight of ");
		// printf("%f\n", inactive_w_lambda[i]);*/
		// group_norms[i] = sqrt(group_norms[i])/inactive_w_lambda[i];
		// /*compare_val = lambda[0]*inactive_w_lambda[i];/*
		// printf("lambda val:");
		// printf("%f\n",lambda[0]);
		// printf("lambda weight:");
		// printf("%f\n",inactive_w_lambda[i]);
		// printf("comapre val:");
		// printf("%f\n", compare_val);
		// printf("i val is: ");
		// printf("%d\n",i);
		// printf("norm value is ");
		// printf("%f\n", sqrt(group_norms[i]));*/
		// if(group_norms[i] > lambda[0]){
			// //printf("working on the failure problem");
			// fail_idx[i] = 1;
			// //printf("fail idx value is:");
			// //printf("%d\n",fail_idx[i]);
			// *kkt_fail = 1;
			// //printf("kkt_fail value is:");
			// //printf("%d\n", kkt_fail[0]);
			// //printf("done working on the failure problem");
		// } else{
			// //printf("not a failure");
			// fail_idx = 0;
		// }
	// }
// }





/*
struct groupMat{
	int size;	
	double vals[];
};*/

// void jaggedArray(double* x, int* groups, int* n, int* p, int* g_num, int* g_count){
	// int i, j, k, m; 
	// struct groupMat ja[g_num[0]];
	// //double recentArray[16*40];
	// for(i = 0; i < g_num[0]; i++){
		// //printf("Working on group %d", i);
		// //printf(" With entries of: %d \n", g_count[i]);
		// int groupSize = g_count[i]*n[0];
		// double temp_array[groupSize];
		// int array_pos = 0; 
		// for(j = 0; j < p[0]; j++){
			// //printf("Checking predictor %d \n", j);
			// if(groups[j]==(i+1)){
				// //printf("Predictor in Group %d %d \n", j, i);
				// for(k = 0; k < n[0]; k++){
					// //printf("X entry is \%d \n", j*n[0]+k);
					// temp_array[array_pos] = x[j*n[0]+k];
					// //printf("Group Position Entry %d %d %f \n", i, array_pos, temp_array[array_pos]);
					// array_pos++; 
				// }
			// }
		// }
		// //printf("attempting to store array for group %d \n", i);
		// //ja[i].vals = temp_array;
		// ja[i].size = groupSize;
		// /*printf("looping through temp array \n");
		// for(m=0; m < groupSize; m++){
			// printf("Group Entry Value %d %d %f \n", i, m, temp_array[m]);
		// }*/
		// if(i==3){
			// printf("looping through vals \n");
			// for(m=0; m < groupSize; m++){
				// printf("Group Entry Value %d %d %f \n", i, m, ja[i].vals[m]);
			// }
		// }
		// //recentArray = temp_array;
	// }
	// i=3;
	// printf("looping through vals \n");
	// for(m=0; m < ja[i].size; m++){
		// printf("Outside Loop Group Entry Value %d %d %f \n", i, m, ja[i].vals[m]);
	// }
	
	// /*for(i = 0; i < g_num[0]; i++){
		// for(m=0; m < ja[i].size; m++){
				// printf("Second Time Group Entry Value %d %d %f \n", i, m, ja[i].vals[m]);
		// }
	// }
	// for(i = 0; i < g_num[0]; i++){
		// for(m=0; m < ja[i].size; m++){
				// printf("Third Time Group Entry Value %d %d %f \n", i, m, ja[i].vals[m]);
		// }
	// }*/
	// /*
	// //Now we will try to loop through the jagged Array
	// printf("Looping through the jagged array again \n");
	// for(i = 0; i < g_num[0]; i++){
		// int gs = ja[i].size; 
		// printf("Group Size %d %d \n", i, gs);
		// for(j=0; j < gs; j++){
			// double val = ja[i].vals[j]; 
			// printf("Group Entry Value %d %d %f \n", i, j, val);
		// }
	// }*/
	
	// i = 0;
	// for(m=0; m < ja[i].size; m++){
			// printf("Third Time Group Entry Value %d %d %f \n", i, m, ja[i].vals[m]);
	// }
	
	// i = 3;
	// for(m=0; m < ja[i].size; m++){
			// printf("Third Time Group Entry Value %d %d %f \n", i, m, ja[i].vals[m]);
	// }*/
// }

void getGroupCol(double* x, int* n, int* cols, int col_num, double* output){
	int i, j, k;
	k = 0; 
	//printf("column num is %d \n", col_num);
	for(i = 0; i < col_num; i++){
		//printf("Working on colum %d \n", cols[i]);
		for(j = 0; j < n[0]; j++){
			//printf("Inputing entry %d \n", cols[i]*n[0]+j);
			//printf("With values of %d %d %d \n", cols[i],n[0],j);
			output[k] = x[cols[i]*n[0]+j];
			k++;
		}
	}
}



void update_residuals(double* r, double* x, int col, int* n, double beta_new, double* beta_old){
	int i; 
	for(i = 0; i < n[0]; i++){
		r[i] = r[i] - x[col*n[0]+i] * (beta_new-beta_old[col]);
	}
}


//residuals need to be provided. 
void solve_beta(double* y, double* x, double* tau, double* gamma, 
			    double* weights, double* lambdaj, 
				double* w_lambda, double* eigenval, double* beta, int* max_iter, double* epsilon, int* apprx, int* n, int* p, int* ng,
				int* gn, int* converge, double* r, int* iter){
  
	int k;
	double delta, new_val;
	double u_norm[1];
	iter[0] = 0;
	delta = 2;
	double neg_grad[1]; 
	//int update_track;
	//update_track = 0;
  while (delta > epsilon[0] && iter[0]<max_iter[0]){   
    iter[0] = iter[0]+1;
	//printf("In loop for lambda iter %f %d \n", lambdaj[0], iter[0]);
	neg_grad[0] = 0; 
/*	print_dbl_array(r, n[0]);
	print_dbl_array(weights, n[0]);
	print_dbl_array(tau, );
	print_dbl_array(gamma);
	print_dbl_array(x);
	print_array(n);
	print_array(apprx);
	print_dbl_array(neg_grad);*/
	neg_gradient(r, weights, tau, gamma, x, n, 1, apprx, 0, neg_grad);
	//printf("Negative gradient value is %f \n", neg_grad[0]);
	//printf("Old Intercept is %f \n", beta[0]);
	//if(apprx[0]==1){
	//printf("Old intercept is %f \n", beta[0]);
	new_val = beta[0]+neg_grad[0]*gamma[0];
	//printf("Intercept adjustment is %f \n", neg_grad[0]*gamma[0]);
	//} else{
	//	beta[0] = beta[0]+neg_grad[0];
	//}
	//printf("New intercept is %f \n", beta[0]);
	delta = fabs(neg_grad[0]);

    //vec_scalar_subtract(r,neg_grad,n);
	update_residuals(r, x, 0, n, new_val, beta);
	beta[0] = new_val;
	
	//printf("residuals after intercept update \n");
	//print_dbl_array(r,n[0]);
	//printf("residual one is %f \n", r[0]);

    
    // update for each group
	//int start_pos = n[0]; 
    int col_pos = 1; // 0 is the intercept
	for(k = 0; k < ng[0]; k++){		
	//	printf("working on group %d \n", k);
		int group_n = gn[k];
	//	printf("initializing cols with %d \n", gn[k]);
	/*
		int cols[group_n];
		int i, col_spot; 
		for(i = 0; i < p[0]; i++){
			if(group_index[i]== k+1){
				cols[col_spot] = i;
				//printf("cols pos and value %d %d %d %d \n", col_spot, i, cols[col_spot], cols[col_spot]+1);
			}
		}*/
		// This ends the code for the above comment
		double u[group_n];
		u_norm[0] = 0;

		
		//printf("neg gradient for group %d \n", k+1);
		/*printf("r one %f \n", r[0]);
		printf("weights one %f \n", weights[0]);
		printf("tau %f \n", tau[0]);
		printf("gamma %f \n", gamma[0]);
		printf("x start %f \n", x[start_pos]);
		printf("n value %d \n", n[0]);
		printf("num in group %d \n", gn[k]);
		printf("apprx value is %d \n", apprx[0]);
		printf("start pos is %d \n", start_pos);*/
		neg_gradient(r, weights, tau, gamma, x, n, gn[k], apprx, col_pos*n[0], u);
		//start_pos = start_pos+n[0]*gn[k];
		
		//printf("neg gradient calculated \n");
		//printf("cols pos and value %d %d \n", 0, cols[0]);
		/*int m; 
		for(m =0; m < gn[k]; m++){
			printf("neg gradient pos value %d %f \n", m, u[m]);
		}*/
		int i; 
		//printf("about to update u for this many betas: %d \n", gn[k]);
		for(i = 0; i < gn[k]; i++){
			//printf("cols pos and value %d %d %d \n", i, cols[i], cols[i]+1);
			//printf("working on index betaEntry %d %d \n", i, cols[i]+1);
			//printf("adding together u eigenvalue beta %f %f %f \n", u[i], eigenval[k], beta[col_pos+i]);
			u[i] = u[i] + eigenval[k]*beta[col_pos+i];
		}
		//printf("u value below for group of size current_u_norm %d %f \n", gn[k], u_norm[0]);
		//print_dbl_array(u,group_n);
		l2norm(u, gn[k], u_norm);
		//printf("paste u_norm and compare %f %f \n", u_norm[0], lambdaj[0]*w_lambda[k]);
		//printf("about to update residuals \n");
		if(u_norm[0] <= lambdaj[0]*w_lambda[k]){
			for(i = 0; i < gn[k]; i++){
				//printf("old beta pos value %d %f \n", col_pos+i, beta[col_pos+i]);
				update_residuals(r, x, col_pos+i, n, 0, beta);
			//	print_dbl_array(r,n[0]);
				if(fabs(beta[col_pos+i]) > delta){
					delta = fabs(beta[col_pos+i]);
				}
				beta[col_pos+i] = 0; 		
				//printf("new beta pos value %d %f \n", col_pos+i, beta[col_pos+i]);				
			}
		} else{
			for(i = 0; i < gn[k]; i++){
				//printf("beta update values temp1 lambda w_lambda u_norm eigenval %f %f %f %f %f \n", u[i], lambdaj[0], w_lambda[k], u_norm[0], eigenval[k]);
				//printf("old beta pos value %d %f \n", col_pos+i, beta[col_pos+i]);
				new_val = u[i]*(1-lambdaj[0]*w_lambda[k]/u_norm[0])/eigenval[k];
				update_residuals(r, x, col_pos+i, n, new_val, beta);
		//		print_dbl_array(r,n[0]);
				if(fabs(beta[col_pos+i]-new_val) > delta){
					delta = fabs(beta[col_pos+i]-new_val);
				}
				beta[col_pos+i] = new_val;
				//printf("new beta pos value %d %f \n", col_pos+i, beta[col_pos+i]);				
			}
		}
		col_pos = col_pos + gn[k];
		// stopped here. 
		//#r0<- as.vector(y-int1-as.matrix(x)%*%beta1)
		//r0<- as.vector(r0-as.matrix(x[,ind])%*%(beta1_k-beta0[ind]))      
	//	printf("residuals after intercept update to group %d \n", k+1);
	//	print_dbl_array(r,n[0]);

    }
    
    //delta<- max(abs(c(int1, beta1)- c(int0, beta0)))
    //#delta.seq<- c(delta.seq, delta)
   // printf("iter value delta value %d %f \n", iter[0], delta);
  }
  
  if(iter[0] < max_iter[0]){
	  converge[0] = 1;
  } else{
	  converge[0] = 0; 
  }
					
}

