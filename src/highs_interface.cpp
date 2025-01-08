
#include "highs_interface.h"


using namespace Rcpp;


/*
class HighsR: public Highs {
    public:

    void set_ncol(int32_t ncol) {
        model_.lp_.num_col_ = ncol;
    }

};
*/


static void R_message_handler(HighsLogType type, const char* message, void* log_callback_data) {
  Rcpp::Rcout << message << std::endl;
}


std::vector<HighsVarType> to_vartype(std::vector<int32_t> type) {
    std::vector<HighsVarType> htype;
    std::vector<HighsVarType> variable_types = {
        HighsVarType::kContinuous, HighsVarType::kInteger, HighsVarType::kSemiContinuous,
        HighsVarType::kSemiInteger, HighsVarType::kImplicitInteger};
    for (std::size_t i = 0; i < type.size(); ++i) {
        htype.push_back(variable_types[type[i]]);
    }
    return htype;
}


// [[Rcpp::export]]
SEXP new_model() {
    Rcpp::XPtr<HighsModel> highs_model(new HighsModel(), true);
    return highs_model;
}

// [[Rcpp::export]]
SEXP model_set_ncol(SEXP mpt, int32_t ncol) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.num_col_ = ncol;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_nrow(SEXP mpt, int32_t nrow) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.num_row_ = nrow;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_sense(SEXP mpt, bool maximum) {
    Rcpp::XPtr<HighsModel>model(mpt);
    if (maximum) {
        model->lp_.sense_ = ObjSense::kMaximize;
    } else {
        model->lp_.sense_ = ObjSense::kMinimize;
    }
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_offset(SEXP mpt, double_t offset) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.offset_ = offset;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_objective(SEXP mpt, std::vector<double> objective) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.col_cost_ = objective;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_lower(SEXP mpt, std::vector<double> lower) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.col_lower_ = lower;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_upper(SEXP mpt, std::vector<double> upper) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.col_upper_ = upper;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_constraint_matrix_(SEXP mpt, std::string format,
                                  std::vector<int32_t> start,
                                  std::vector<int32_t> index,
                                  std::vector<double> value) {
    Rcpp::XPtr<HighsModel>model(mpt);
    if (format == "colwise") {
        model->lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    } else if (format == "rowwise") {
        model->lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
    } else if (format == "rowwise_partitioned") {
        model->lp_.a_matrix_.format_ = MatrixFormat::kRowwisePartitioned;
    } else {
        Rcpp::stop("unkown format!");
    }
    model->lp_.a_matrix_.start_ = start;
    model->lp_.a_matrix_.index_ = index;
    model->lp_.a_matrix_.value_ = value;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_lhs(SEXP mpt, std::vector<double> lower) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.row_lower_ = lower;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_rhs(SEXP mpt, std::vector<double> upper) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->lp_.row_upper_ = upper;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_hessian_(SEXP mpt, std::string format, int32_t dim,
                        std::vector<int32_t> start,
                        std::vector<int32_t> index,
                        std::vector<double> value) {
    Rcpp::XPtr<HighsModel>model(mpt);
    model->hessian_.dim_ = dim;
    if (format == "triangular") {
        model->hessian_.format_ = HessianFormat::kTriangular;
    } else if (format == "square") {
        model->hessian_.format_ = HessianFormat::kSquare;
    } else {
        Rcpp::stop("unkown format!");
    }
    model->hessian_.start_ = start;
    model->hessian_.index_ = index;
    model->hessian_.value_ = value;
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP model_set_vartype(SEXP mpt, std::vector<int32_t> type) {
    Rcpp::XPtr<HighsModel>model(mpt);
    if (model->lp_.integrality_.size() < type.size()) {
        model->lp_.integrality_.resize(type.size());
    }
    std::vector<HighsVarType> variable_types = {
        HighsVarType::kContinuous, HighsVarType::kInteger, HighsVarType::kSemiContinuous,
        HighsVarType::kSemiInteger, HighsVarType::kImplicitInteger};
    for (std::size_t i = 0; i < type.size(); ++i) {
        model->lp_.integrality_[i] = variable_types[type[i]];
    }
    return R_NilValue;
}

// [[Rcpp::export]]
int32_t model_get_nvars(SEXP mpt) {
    Rcpp::XPtr<HighsModel>model(mpt);
    return static_cast<int32_t>(model->lp_.num_col_);
}

// [[Rcpp::export]]
int32_t model_get_ncons(SEXP mpt) {
    Rcpp::XPtr<HighsModel>model(mpt);
    return static_cast<int32_t>(model->lp_.num_row_);
}


RCPP_MODULE(RcppHighs) {
    class_<Highs>("Highs")

    .constructor()

    
    .method("getObjectiveValue", &Highs::getObjectiveValue)
    .method("getNumCol", &Highs::getNumCol)
    .method("getNumRow", &Highs::getNumRow)
    .method("getNumNz", &Highs::getNumNz)
    .method("getHessianNumNz", &Highs::getHessianNumNz)
    // .method("getObjectiveSense", &Highs::getObjectiveSense)
    // .method("getObjectiveOffset", &Highs::getObjectiveOffset)

    ;
}


// [[Rcpp::export]]
SEXP new_solver(SEXP mpt) {
    Rcpp::XPtr<Highs> highs(new Highs(), true);
    highs->setLogCallback(R_message_handler);
    if (Rf_isNull(mpt)) {
        return highs;
    }
    Rcpp::XPtr<HighsModel>model(mpt);
    HighsStatus return_status = highs->passModel(*model.get());
    if (return_status != HighsStatus::kOk) {
        return R_NilValue;
    }
    return highs;
}


// [[Rcpp::export]]
SEXP highs_pass_model(SEXP hi,
                   const int32_t num_col, const int32_t num_row,
                   const int32_t num_nz, const int32_t a_format,
                   const int32_t sense, const double_t offset,
                   NumericVector col_cost, NumericVector col_lower,
                   NumericVector col_upper, NumericVector row_lower,
                   NumericVector row_upper, IntegerVector a_start,
                   IntegerVector a_index, NumericVector a_value,
                   IntegerVector integrality) {
    Rcpp::XPtr<Highs>highs(hi);
    highs->passModel(num_col,
                     num_row,
                     num_nz,
                     a_format,
                     sense,
                     offset,
                     &(col_cost[0]),
                     &(col_lower[0]),
                     &(col_upper[0]),
                     &(row_lower[0]),
                     &(row_upper[0]),
                     &(a_start[0]),
                     &(a_index[0]),
                     &(a_value[0]),
                     &(integrality[0]));
    return R_NilValue;
}


// [[Rcpp::export]]
SEXP solver_pass_hessian() {
    return R_NilValue;
}


// [[Rcpp::export]]
SEXP solver_pass_constraints () {
    return R_NilValue;
}


// [[Rcpp::export]]
int32_t solver_get_sense(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    ObjSense sense;
    HighsStatus status = highs->getObjectiveSense(sense);
    if (status == HighsStatus::kOk) {
        return sense == ObjSense::kMaximize;
    } else {
        Rcpp::stop("could not obtain the sense of the objective function.");
    }
}

// [[Rcpp::export]]
int32_t solver_set_sense(SEXP hi, bool maximum) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status;
    if (maximum) {
        status = highs->changeObjectiveSense(ObjSense::kMaximize);
    } else {
        status = highs->changeObjectiveSense(ObjSense::kMinimize);
    }
    return static_cast<int32_t>(status);
}

// [[Rcpp::export]]
int32_t solver_set_offset(SEXP hi, const double ext_offset) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->changeObjectiveOffset(ext_offset);
    return static_cast<int32_t>(status);
}

// [[Rcpp::export]]
int32_t solver_set_integrality(SEXP hi, std::vector<int32_t> index, std::vector<int32_t> type) {
    Rcpp::XPtr<Highs>highs(hi);
    const int32_t num_set_entries = index.size();
    const HighsInt* set = &(index[0]);
    const HighsVarType* integrality = &(to_vartype(type)[0]);
    HighsStatus status = highs->changeColsIntegrality(num_set_entries, set, integrality);
    return static_cast<int32_t>(status);
}

// [[Rcpp::export]]
int32_t solver_set_objective(SEXP hi, std::vector<int32_t> index, std::vector<double_t> obj) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->changeColsCost(index.size(), &(index[0]), &(obj[0]));
    return static_cast<int32_t>(status);
}


// [[Rcpp::export]]
int32_t solver_set_variable_bounds(SEXP hi, std::vector<int32_t> index, std::vector<double_t> lower, std::vector<double_t> upper) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->changeColsBounds(index.size(), &(index[0]), &(lower[0]), &(upper[0]));
    return static_cast<int32_t>(status);
}


// [[Rcpp::export]]
int32_t solver_set_constraint_bounds(SEXP hi, std::vector<int32_t> index, std::vector<double_t> lower, std::vector<double_t> upper) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->changeRowsBounds(index.size(), &(index[0]), &(lower[0]), &(upper[0]));
    return static_cast<int32_t>(status);
}


// [[Rcpp::export]]
SEXP solver_set_coeff(SEXP hi, std::vector<int32_t> row, std::vector<int32_t> col, std::vector<double_t> val) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status;
    for (std::size_t i = 0; i < row.size(); ++i) {
        status = highs->changeCoeff(row[i], col[i], val[i]);
        if (status != HighsStatus::kOk) {
            Rcpp::stop("could not change the coefficient.");
        }
    }
    return R_NilValue;
}


// [[Rcpp::export]]
int32_t solver_add_vars(SEXP hi, std::vector<double_t> lower, std::vector<double_t> upper) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status = highs->addVars(lower.size(), &(lower[0]), &(upper[0]));
    return static_cast<int32_t>(status);
}

/*
model_set_lower
model_set_upper
model_set_constraint_matrix_
model_set_lhs
model_set_rhs
model_set_hessian_
model_set_vartype
*/
//
// changeObjectiveSenseInterface
// 
// addRowsInterface
// addColsInterface
// changeColsIntegrality
// changeColsCost
// changeColsBounds
// changeRowsBounds
// changeCoeff
// addCols
// addRows
// setSolution
// getHotStart
// setHotStart

// [[Rcpp::export]]
int32_t solver_set_option(SEXP hi, std::string key, SEXP value) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus status;
    if (Rf_isLogical(value)) {
        bool logval = Rcpp::as<bool>(value);
        status = highs->setOptionValue(key, logval);
    } else if (Rf_isInteger(value)) {
        HighsInt intval = Rcpp::as<int32_t>(value);
        status = highs->setOptionValue(key, intval);
    } else if (Rf_isNumeric(value)) {
        double numval = Rcpp::as<double>(value);
        status = highs->setOptionValue(key, numval);
    } else if (Rf_isString(value)) {
        std::string strval = Rcpp::as<std::string>(value);
        status = highs->setOptionValue(key, strval);
    } else {
        Rcpp::stop("unkown type of value.");
    }
    return static_cast<int32_t>(status);
}


// [[Rcpp::export]]
int32_t solver_clear(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->clear();
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_clear_model(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->clearModel();
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_clear_solver(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->clearSolver();
    return static_cast<int32_t>(return_status);
}



// [[Rcpp::export]]
int32_t solver_run(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->run();
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
SEXP solver_get_model(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsModel mpt = highs->getModel();
    Rcpp::XPtr<HighsModel> highs_model(&mpt, true);
    return highs_model;
}


// [[Rcpp::export]]
int32_t solver_get_num_col(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    return highs->getNumCol();
}


// [[Rcpp::export]]
int32_t solver_get_num_row(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    return highs->getNumRow();
}



//
// int32_t solver_presolve(SEXP hi) {
//     Rcpp::XPtr<Highs>highs(hi);
//     HighsStatus return_status = highs->presolve();
//     return static_cast<int32_t>(return_status);
// }


// TODO
// HighsStatus postsolve(const HighsSolution& solution, const HighsBasis& basis);


// TODO
// const HighsLp& getPresolvedLp()


// TODO
// const HighsModel& getPresolvedModel()


// [[Rcpp::export]]
int32_t solver_write_model(SEXP hi, const std::string filename) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->writeModel(filename);
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_write_basis(SEXP hi, const std::string filename) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->writeBasis(filename);
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
std::string solver_status_message(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsModelStatus& model_status = highs->getModelStatus();
    return highs->modelStatusToString(model_status);
}


// [[Rcpp::export]]
int32_t solver_status(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsModelStatus& model_status = highs->getModelStatus();
    return static_cast<int32_t>(model_status);
}


// [[Rcpp::export]]
double_t  solver_infinity() {
    Highs highs;
    double_t dbl_inf = highs.getInfinity();
    return dbl_inf;
}


// [[Rcpp::export]]
SEXP reset_global_scheduler(bool blocking) {
    Highs highs;
    highs.resetGlobalScheduler(blocking);
    return R_NilValue;
}


// [[Rcpp::export]]
Rcpp::List solver_info(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsInfo& info = highs->getInfo();
    List z = List::create(
        Named("valid") = info.valid,
        Named("mip_node_count") = info.mip_node_count,
        Named("simplex_iteration_count") = info.simplex_iteration_count,
        Named("ipm_iteration_count") = info.ipm_iteration_count,
        Named("qp_iteration_count") = info.qp_iteration_count,
        Named("crossover_iteration_count") = info.crossover_iteration_count,
        Named("primal_solution_status") = highs->solutionStatusToString(info.primal_solution_status),
        Named("dual_solution_status") = highs->solutionStatusToString(info.dual_solution_status),
        Named("basis_validity") = info.basis_validity,
        Named("objective_function_value") = info.objective_function_value,
        Named("mip_dual_bound") = info.mip_dual_bound,
        Named("mip_gap") = info.mip_gap,
        Named("num_primal_infeasibilities") = info.num_primal_infeasibilities,
        Named("max_primal_infeasibility") = info.max_primal_infeasibility,
        Named("sum_primal_infeasibilities") = info.sum_primal_infeasibilities,
        Named("num_dual_infeasibilities") = info.num_dual_infeasibilities,
        Named("max_dual_infeasibility") = info.max_dual_infeasibility,
        Named("sum_dual_infeasibilities") = info.sum_dual_infeasibilities
    );
    return z;
}

// [[Rcpp::export]]
Rcpp::List solver_solution(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsSolution& solution = highs->getSolution();
    List z = List::create(
        Named("value_valid") = solution.value_valid,
        Named("dual_valid") = solution.dual_valid,
        Named("col_value") = solution.col_value,
        Named("col_dual") = solution.col_dual,
        Named("row_value") = solution.row_value,
        Named("row_dual") = solution.row_dual
    );
    return z;
}

// [[Rcpp::export]]
bool solver_get_bool_option(SEXP hi, std::string key) {
    Rcpp::XPtr<Highs>highs(hi);
    bool value;
    highs->getOptionValue(key, value);
    // HighsInt    
    return value;
}

// [[Rcpp::export]]
int32_t solver_get_int_option(SEXP hi, std::string key) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsInt value;
    highs->getOptionValue(key, value);
    return value;
}

// [[Rcpp::export]]
double_t solver_get_dbl_option(SEXP hi, std::string key) {
    Rcpp::XPtr<Highs>highs(hi);
    double_t value;
    highs->getOptionValue(key, value);
    return value;
}

// [[Rcpp::export]]
std::string solver_get_str_option(SEXP hi, std::string key) {
    Rcpp::XPtr<Highs>highs(hi);
    std::string value;
    highs->getOptionValue(key, value);   
    return value;
}


// [[Rcpp::export]]
int32_t solver_change_variable_bounds(SEXP hi, IntegerVector idx, NumericVector lower, NumericVector upper) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->changeColsBounds(idx.size(), &(idx[0]), &(lower[0]), &(upper[0]));
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_change_constraint_bounds(SEXP hi, IntegerVector idx, NumericVector lhs, NumericVector rhs) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->changeRowsBounds(idx.size(), &(idx[0]), &(lhs[0]), &(rhs[0]));
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_add_rows(SEXP hi, NumericVector lhs, NumericVector rhs,
    IntegerVector start, IntegerVector index, NumericVector value) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->addRows(
        lhs.size(), //!< num_new_row = Number of new rows
        &(lhs[0]),  //!< lower = Array of size num_new_row with lower bounds
        &(rhs[0]),  //!< upper = Array of size num_new_row with upper bounds
        value.size(), //!< num_new_nz = Number of new nonzeros
        &(start[0]),   //!< starts = Array of size num_new_row with start indices of the rows
        &(index[0]),   //!< indices = Array of size num_new_nz with column indices for all rows
        &(value[0])    //!< values = Array of size num_new_nz with column values for all rows
    );
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
int32_t solver_add_cols(SEXP hi, NumericVector costs,
    NumericVector lower, NumericVector upper,
    IntegerVector start, IntegerVector index, NumericVector value) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsStatus return_status = highs->addCols(
        lower.size(), //!< num_new_col = Number of new columns
        &(costs[0]),   //!< costs = Array of size num_new_col with costs
        &(lower[0]),   //!< lower = Array of size num_new_col with lower bounds
        &(upper[0]),   //!< upper = Array of size num_new_col with upper bounds
        value.size(), //!< num_new_nz = Number of new nonzeros
        &(start[0]),   //!< starts = Array of size num_new_row with start indices of the columns
        &(index[0]),   //!< indices = Array of size num_new_nz with row indices for all columns
        &(value[0])    //!< values = Array of size num_new_nz with row values for all columns
    );
    return static_cast<int32_t>(return_status);
}


// [[Rcpp::export]]
Rcpp::NumericVector solver_get_lp_costs(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    int32_t nvar = highs->getNumCol();
    NumericVector costs(nvar);
    HighsModel model = highs->getModel();
    for (int32_t i = 0; i < nvar; i++) {
        costs[i] = model.lp_.col_cost_[i];
    }
    return costs;
}


// [[Rcpp::export]]
Rcpp::NumericVector solver_get_variable_bounds(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    int32_t nvar = highs->getNumCol();
    NumericVector vbounds(2 * nvar);
    HighsModel model = highs->getModel();
    for (int32_t i = 0; i < nvar; i++) {
        vbounds[i] = model.lp_.col_lower_[i];
        vbounds[nvar + i] = model.lp_.col_upper_[i];
    }
    return vbounds;
}


// [[Rcpp::export]]
Rcpp::NumericVector solver_get_constraint_bounds(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    int32_t nvar = highs->getNumCol();
    NumericVector lhs_rhs(2 * nvar);
    HighsModel model = highs->getModel();
    for (int32_t i = 0; i < nvar; i++) {
        lhs_rhs[i] = model.lp_.row_lower_[i];
        lhs_rhs[nvar + i] = model.lp_.row_upper_[i];
    }
    return lhs_rhs;
}


// [[Rcpp::export]]
Rcpp::List solver_get_constraint_matrix(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    HighsModel model = highs->getModel();
    HighsSparseMatrix A = model.lp_.a_matrix_;
    List z = List::create(
        Named("format") = static_cast<int32_t>(A.format_),
        Named("nrow") = A.num_row_,
        Named("ncol") = A.num_col_,
        Named("start") = A.start_,
        Named("p_end") = A.p_end_,
        Named("index") = A.index_,
        Named("value") = A.value_
    );
    return z;
}


// [[Rcpp::export]]
Rcpp::IntegerVector solver_get_vartype(SEXP hi) {
    Rcpp::XPtr<Highs>highs(hi);
    const HighsLp& lp = highs->getLp();
    int32_t m = lp.integrality_.size();
    IntegerVector type(m);
    for (int32_t i = 0; i < type.size(); ++i) {
        type[i] = static_cast<int32_t>(lp.integrality_[i]);
    }
    return type;
}
