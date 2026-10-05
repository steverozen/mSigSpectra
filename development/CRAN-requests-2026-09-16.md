Please always explain all acronyms in the description text.
For more details: <https://contributor.r-project.org/cran-cookbook/description_issues.html#explaining-acronyms>

-> SBS, DBS, etc.

Please always write package names, software names and API (application programming interface) names in single quotes in title and description.
e.g: --> 'shiny'
Please note that package names are case sensitive.
For more details: <https://contributor.r-project.org/cran-cookbook/description_issues.html#formatting-software-names>

License components with restrictions and base license permitting such:
GPL-3 + file LICENSE

We do not need "+ file LICENSE" and the file as these are part of R. This is only needed in case of attribution requirements or other possible restrictions. Hence please omit it.
For more details: <https://contributor.r-project.org/cran-cookbook/description_issues.html#license-files>

Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar)
For more details: <https://contributor.r-project.org/cran-cookbook/docs_issues.html#missing-value-tags-in-.rd-files>

-> Missing Rd-tags: 
        check_and_remove_discarded_variants.Rd: \value 
        is_catalog.Rd: \value 
        subset_catalog.Rd: \value


You write information messages to the console that cannot be easily suppressed.
It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object. Instead of print()/cat() rather use message()/warning() or if(verbose)cat(..) (or maybe stop()) if you really have to write text to the console. (except for print, summary, interactive functions)
For more details: <https://contributor.r-project.org/cran-cookbook/code_issues.html#using-printcat>


-> R/quick_check_vcf.R

Please fix and resubmit. Best, Leonore Hochhauser (she/her)

This from an email fom Stored with zero-access encryption
Leonore Hochhauser<leonore.hochhauser@wu.ac.at> dated 2026 09 15

