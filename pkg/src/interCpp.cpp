
#include <Rcpp.h>

RcppExport SEXP interCpp(SEXP set1R, SEXP set2R) {

BEGIN_RCPP

  Rcpp::NumericVector set1(set1R);
  Rcpp::NumericVector set2(set2R);

  std::sort(set1.begin(), set1.end());
  std::sort(set2.begin(), set2.end());

  typedef std::set<int, std::less<int> > intSet;
  intSet s1, s2, s3;
  s1.insert(set1.begin(), set1.end());
  s2.insert(set2.begin(), set2.end());
 
  std::set_intersection( set1.begin(), set1.end(), set2.begin(), set2.end(),
   std::insert_iterator<intSet>(s3,s3.begin()) );

  return Rcpp::List::create(Rcpp::Named("inter") = s3);
  
END_RCPP
}
