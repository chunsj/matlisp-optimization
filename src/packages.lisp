(in-package #:cl-user)

(defpackage "SR1"
  (:nicknames :sr1)
  (:import-from :λ-reader #:λ)
  (:use :cl :matlisp :matlisp-utilities #:iterate)
  (:export #:update! #:l-update! #:l-query #:lsr1 #:sr1-warning))

(defpackage "BFGS"
  (:nicknames #:bfgs)
  (:import-from #:λ-reader #:λ)
  (:use :cl :matlisp :matlisp-utilities #:iterate)
  (:export #:update! #:l-update! #:l-query #:lbfgs #:bfgs-warning))

(defpackage "MATLISP-OPTIMIZATION"
  (:nicknames #:opt)
  (:import-from #:λ-reader #:λ)
  (:use :cl :matlisp :matlisp-utilities #:iterate)
  (:export #:backtracking-linesearch))

(defpackage "MATLISP-OPTIMIZATION/TESTS"
  (:import-from #:λ-reader #:λ)
  (:use :cl :matlisp :matlisp-utilities #:iterate)
  (:export #:matlisp-optimization #:run-tests))

(fiveam:def-suite matlisp-optimization/tests:matlisp-optimization
    :description "Regression tests for matlisp-optimization")
(defun matlisp-optimization/tests:run-tests ()
  (5am:run! 'matlisp-optimization/tests:matlisp-optimization))
(fiveam:in-suite matlisp-optimization/tests:matlisp-optimization)
