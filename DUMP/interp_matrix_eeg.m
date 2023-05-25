function struct_out = interp_matrix_eeg(struct_in,interp_fact_x,interp_fact_y)


X = struct_in.tf_time;
Y = struct_in.tf_frex;
[Xq,Yq] = meshgrid(min(X):interp_fact_x:max(X),min(Y):interp_fact_y:max(Y));
xq = min(X):interp_fact_x:max(X);
yq = min(Y):interp_fact_y:max(Y);

data_in = struct_in.tf_data_subs;
dims_data = size(data_in);
data_out = zeros(dims_data(1),dims_data(2),length(yq),length(xq));


for i = 1:dims_data(1)
    for j = 1:dims_data(2)
        c_dat = squeeze(data_in(i,j,:,:));
        data_out(i,j,:,:) = interp2(X,Y,c_dat,Xq,Yq);
    end
end

struct_out = struct_in ;
struct_out.tf_data_subs = data_out;
struct_out.tf_time = xq;
struct_out.tf_frex = yq;
struct_out.tf_data = squeeze(mean(data_out));


end