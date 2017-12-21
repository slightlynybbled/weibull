import pytest
import weibull


def test_version():
    """Test for the presence of the version as a string."""
    assert '.' in weibull.__version__


def test_non_fitted():
    """
    Tests for proper fail condition if beta and eta do not exist

    :return: None
    """
    analysis = weibull.Analysis([0], [False])

    assert analysis.beta is None
    assert analysis.eta is None

    with pytest.raises(weibull.ParameterError):
        analysis.probplot()

    with pytest.raises(weibull.ParameterError):
        analysis.pdf()

    with pytest.raises(weibull.ParameterError):
        analysis.sf()

    with pytest.raises(weibull.ParameterError):
        analysis.hazard()

    with pytest.raises(weibull.ParameterError):
        analysis.cdf()

    with pytest.raises(weibull.ParameterError):
        analysis.fr()

    with pytest.raises(weibull.ParameterError):
        analysis.b(10)

    with pytest.raises(weibull.ParameterError):
        m = analysis.mean

    with pytest.raises(weibull.ParameterError):
        m = analysis.mttf

    with pytest.raises(weibull.ParameterError):
        m = analysis.median

    with pytest.raises(weibull.ParameterError):
        m = analysis.characteristic_life


def test_complete_dataset():
    """
    Tests basic calculation of beta and eta for a complete data set
    along with some misc. functionality based on the presence of
    beta and eta.

    :return:
    """
    fail_times = [
        9402.7,
        6082.4,
        13367.2,
        10644.6,
        8632.0,
        3043.4,
        12860.2,
        1034.5,
        2550.9,
        3637.1
    ]

    analysis = weibull.Analysis(fail_times)
    analysis.fit()

    assert 1.345 < analysis.beta < 1.355
    assert 8100.0 < analysis.eta < 8160.0
    assert 1500.0 < analysis.b(10) < 1600
    assert 7400.0 < analysis.mttf < 7500
    assert 8100.0 < analysis.eta < 8160.0


def test_censored_dataset():
    current_run_time = 4200.0

    fail_times = [current_run_time] * 10
    fail_times[7] = 1034.5
    fail_times[8] = 2550.9
    fail_times[6] = 3043.4

    suspended = [True, True, True, True, True,
                 False, False, False, True, True]

    analysis = weibull.Analysis(fail_times, suspended=suspended, unit='hour')
    analysis.fit()

    assert 1.5 < analysis.beta < 1.6
    assert 6500.0 < analysis.eta < 6630


def test_design_bad_reliability():
    with pytest.raises(ValueError):
        designer = weibull.Design(
            target_cycles=10000,
            reliability=0.9999,
            confidence_level=0.90,
            expected_beta=1.5
        )


def test_design_bad_confidence():
    with pytest.raises(ValueError):
        designer = weibull.Design(
            target_cycles=10000,
            reliability=0.9,
            confidence_level=0.0001,
            expected_beta=1.5
        )


def test_design():
    designer = weibull.Design(
        target_cycles=10000,
        reliability=0.9,
        confidence_level=0.90,
        expected_beta=1.5
    )

    assert 21.5 < designer.num_of_units(test_cycles=10000) < 22.0
    assert 10500.0 < designer.num_of_cycles(number_of_units=20) < 10700.0


def test_weibayes_bad_confidence():
    N = 62
    H = 62.0e6

    run_times_desired = [H] * N

    with pytest.raises(ValueError):
        weibayes = weibull.Weibayes(run_times_desired, confidence_level=0.9999, beta=2)


def test_weibayes():
    N = 62
    H = 62.0e6

    run_times_desired = [H] * N

    weibayes = weibull.Weibayes(run_times_desired, confidence_level=0.95, beta=2)
    assert 40050430.0 < weibayes.b(2) < 40100000.0


def test_weibayes_bad_b_params():
    N = 62
    H = 62.0e6

    run_times_desired = [H] * N

    weibayes = weibull.Weibayes(run_times_desired, confidence_level=0.95, beta=2)

    with pytest.raises(ValueError):
        b = weibayes.b(0.1)

    with pytest.raises(ValueError):
        b = weibayes.b(2, 0.999)
