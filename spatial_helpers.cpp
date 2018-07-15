//spatial_helpers

#include<cmath>
#include<algorithm>
#include<iostream>

#define dbg(x) x

struct location{
	double lat, lng;
	location(double a, double b){
        lat=a; lng=b;
	}
	location(){lat=91,lng=181;}
};

const double PI = 3.141592653589793;

inline double radians(double a){return a * PI/180;}
inline double degrees(double a){return a * 180/PI;}

location np_mean(location *locs, int len){
    double a,b;
    a=b=0;
    for(int i=0;i<len;++i){
        a+=locs[i].lat;
        b+=locs[i].lng;
    }
    a/=len;
    b/=len;
    return location(a,b);
}

location compute_midpoint(location pt1, location pt2){
    //convert to radians
    double latitude1, latitude2, latitude3;
    double longitude1, longitude2, longitude3;
    latitude1 = radians(pt1.lat);
    longitude1 = radians(pt1.lng);
    latitude2 = radians(pt2.lat);
    longitude2 = radians(pt2.lng);

    double bx, by;
    bx = cos(latitude2) * cos(longitude2 - longitude1);
    by = cos(latitude2) * sin(longitude2 - longitude1);

    latitude3 = atan2(sin(latitude1) + sin(latitude2), \
           sqrt((cos(latitude1) + bx) * (cos(latitude1) \
           + bx) + by*by));

    return location(degrees(latitude3),degrees(longitude3));
}

double haversine_distance(location point1, location  point2, double radius=6371){
    /*
        (approximate) distance between two points on earth's surface (in KM)
    */
    double latitude1, latitude2, dlatitude;
    double longitude1, longitude2, dlongitude;

    latitude1 = radians(point1.lat);
    longitude1 = radians(point1.lng);
    latitude2 = radians(point2.lat);
    longitude2 = radians(point2.lng);

    dlongitude = longitude2 - longitude1;
    dlatitude = latitude2 - latitude1;

    double a = pow(sin(dlatitude/2),2) + cos(latitude1) * cos(latitude2) * pow(sin(dlongitude/2),2);
    double c = 2 * asin(sqrt(a));
    double km = radius * c;
    return km;
}


location spatial_median(location *X, int len, double eps=1e-3, int max_iter=1000){
    /*
        compute spatial median of a set of > 2 points

        That is, given a set X find the point m s.t.

            sum([dist(x, m) for x in X])

        is minimized. This is a robust estimator of the mode of the set.
    */
    dbg(std::cout<<"inside spatial median\n";)
    int i;
    double *D = new double[len];
    double *Dinv = new double[len];
    double *W = new double[len];
    double r,rinv;
    double Dinvs = 0;
    int iter_ = 0;
    int nonzeros = 0;
    int num_zeros = 0;
    location T,y1,out,R,tmp;
    location y = np_mean(X, len);
    while(true){
        iter_ += 1;
        dbg(std::cout<<iter_<<"\n";)
        Dinvs = 0;
        T.lat=0;
        T.lng=0;
        nonzeros = 0;

        for(i=0; i<len; ++i){
            D[i] = haversine_distance(X[i], y);
            Dinv[i] = D[i]==0 ? 0 : 1/D[i];
            nonzeros = D[i]!=0 ? nonzeros+1 : nonzeros;
            Dinvs += Dinv[i];
        }
        for(i=0; i<len; ++i){
        	W[i] = Dinv[i]/Dinvs;
        	T.lat += W[i] * X[i].lat;
        	T.lng += W[i] * X[i].lng;
        }

		num_zeros = len - nonzeros;

		if(num_zeros == 0)
			y1 = T;
		else if(num_zeros == len)
			return y;
		else{
			R.lat = (T.lat - y.lat) * Dinvs;
			R.lng = (T.lng - y.lng) * Dinvs;
			r = sqrt(R.lat*R.lat + R.lng*R.lng);
			rinv = r==0 ? 0 : num_zeros/r;
			y1.lat = std::max(0.0,1-rinv)*T.lat + std::min(1.0, rinv)*y.lat;
			y1.lng = std::max(0.0,1-rinv)*T.lng + std::min(1.0, rinv)*y.lng;
		}

		tmp.lat = y.lat-y1.lat;
		tmp.lng = y.lng-y1.lng;

        //if (np.linalg.norm(y - y1) < eps) or (iter_ > max_iter)
		if((sqrt(tmp.lat*tmp.lat + tmp.lng*tmp.lng))<eps || (iter_ > max_iter))
			return y1;

		y = y1;

		/*
		D         = np.array(map(lambda x: f(x, y), X))[:,np.newaxis]
        nonzeros  = (D != 0)[:, 0]
        Dinv      = 1 / D[nonzeros]
        Dinvs     = np.sum(Dinv)
        W         = Dinv / Dinvs
        T         = np.sum(W * X[nonzeros], 0)
        num_zeros = len(X) - np.sum(nonzeros)

        if num_zeros == 0:
            y1 = T
        elif num_zeros == len(X):
            out = y
            break
        else:
            R    = (T - y) * Dinvs
            r    = np.linalg.norm(R)
            rinv = 0 if r == 0 else num_zeros/r
            y1   = max(0, 1-rinv)*T + min(1, rinv)*y

        if (np.linalg.norm(y - y1) < eps) or (iter_ > max_iter):
            out = y1
            break

        y = y1*/
    }
}

location spatial_center(location *locs, int len){
    /*
        Computes center of a set of points

        If only one point, "center" is point itself
        If two points, "center" is the midpoint on the line connecting the points
        If > two points, "center" is the spatial median of the set
    */

    if(len==1)
        return locs[0];
    else if(len==2)
        return compute_midpoint(locs[0], locs[1]);
    else
        return spatial_median(locs,len);
}

int main(){

    location l[] = {location(10,10),location(20,20),location(30,30),location(40,40),location(25,25)};
    location ans = spatial_center(l,5);
    dbg(std::cout<<ans.lat<<" "<<ans.lng<<"\n";)
    return 0;
}
