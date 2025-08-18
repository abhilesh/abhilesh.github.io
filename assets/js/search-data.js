// get the ninja-keys element
const ninja = document.querySelector('ninja-keys');

// add the home and posts menu items
ninja.data = [{
    id: "nav-about",
    title: "About",
    section: "Navigation",
    handler: () => {
      window.location.href = "/";
    },
  },{id: "nav-publications",
          title: "Publications",
          description: "More details about my publications can be found on my Google Scholar profile.",
          section: "Navigation",
          handler: () => {
            window.location.href = "/publications/";
          },
        },{id: "nav-blog",
          title: "Blog",
          description: "",
          section: "Navigation",
          handler: () => {
            window.location.href = "/blog/";
          },
        },{id: "nav-cv",
          title: "CV",
          description: "A complete PDF version of my CV is available by clicking the icon on the right.",
          section: "Navigation",
          handler: () => {
            window.location.href = "/cv/";
          },
        },{id: "books-the-making-of-the-atomic-bomb",
          title: 'The Making of the Atomic Bomb',
          description: "",
          section: "Books",handler: () => {
              window.location.href = "/books/the-making-of-the-atomic-bomb/";
            },},{id: "news-joined-university-college-london-s-gee-department-as-a-postdoctoral-research-fellow-to-work-with-prof-max-reuter-and-dr-aida-andres-on-sexual-antagonistic-variation-in-drosophila-melanogaster",
          title: 'Joined University College London’s GEE department as a Postdoctoral Research Fellow to work...',
          description: "",
          section: "News",},{id: "news-published-a-minireview-on-multilevel-selection-on-mtdna-in-current-opinion-in-genetics-amp-amp-development",
          title: 'Published a miniReview on Multilevel Selection on mtDNA in Current Opinion in Genetics...',
          description: "",
          section: "News",},{id: "news-selected-for-participation-in-the-genetics-society-of-america-s-peer-review-training-program-i-will-be-serving-as-a-reviewer-for-the-journals-genetics-and-g3-genes-genomes-genetics",
          title: 'Selected for participation in the Genetics Society of America’s Peer Review Training Program....',
          description: "",
          section: "News",},{id: "news-participated-in-the-embo-population-genomics-background-and-tools-course-in-naples-italy-presented-a-poster-on-the-use-of-approximate-bayesian-computation-in-detecting-sexually-antagonistic-selection-lots-of-in-depth-discussions-on-the-use-of-population-genomic-tools-in-understanding-complex-evolutionary-processes",
          title: 'Participated in the EMBO Population Genomics: Background and Tools Course in Naples, Italy....',
          description: "",
          section: "News",},{id: "news-joined-the-mrc-mitochondrial-biology-unit-at-the-university-of-cambridge-as-a-bioinformatician-i-ll-be-working-with-dr-jelle-van-den-ameele-prof-rita-horvath-and-prof-patrick-chinnery-on-analyzing-mitochondrial-heteroplasmy-and-its-role-in-human-mitochondrial-diseases",
          title: 'Joined the MRC Mitochondrial Biology Unit at the University of Cambridge as a...',
          description: "",
          section: "News",},{
        id: 'social-display_in_header',
        title: 'Display_in_header',
        section: 'Socials',
        handler: () => {
          window.open("", "_blank");
        },
      },{
        id: 'social-email',
        title: 'email',
        section: 'Socials',
        handler: () => {
          window.open("mailto:%61%64%32%33%34%37@%63%61%6D.%61%63.%75%6B", "_blank");
        },
      },{
        id: 'social-scholar',
        title: 'Google Scholar',
        section: 'Socials',
        handler: () => {
          window.open("https://scholar.google.com/citations?user=KmALIboAAAAJ", "_blank");
        },
      },{
        id: 'social-github',
        title: 'GitHub',
        section: 'Socials',
        handler: () => {
          window.open("https://github.com/abhilesh", "_blank");
        },
      },{
        id: 'social-linkedin',
        title: 'LinkedIn',
        section: 'Socials',
        handler: () => {
          window.open("https://www.linkedin.com/in/abhilesh-dhawanjewar", "_blank");
        },
      },{
        id: 'social-blogger',
        title: 'Blogger',
        section: 'Socials',
        handler: () => {
          window.open("https://neuroscience.cam.ac.uk/member/ad2347/", "_blank");
        },
      },{
        id: 'social-dblp',
        title: 'DBLP',
        section: 'Socials',
        handler: () => {
          window.open("https://theealingbeaverproject.com/team/abhilesh-dhawanjewar/", "_blank");
        },
      },{
        id: 'social-instagram',
        title: 'Instagram',
        section: 'Socials',
        handler: () => {
          window.open("https://instagram.com/abhilesh7", "_blank");
        },
      },{
        id: 'social-orcid',
        title: 'ORCID',
        section: 'Socials',
        handler: () => {
          window.open("https://orcid.org/0000-0001-7827-1255", "_blank");
        },
      },{
        id: 'social-researchgate',
        title: 'ResearchGate',
        section: 'Socials',
        handler: () => {
          window.open("https://www.researchgate.net/profile/Abhilesh-Dhawanjewar-2/", "_blank");
        },
      },{
        id: 'social-work',
        title: 'Work',
        section: 'Socials',
        handler: () => {
          window.open("https://www.mrc-mbu.cam.ac.uk/mbu-postdoctoral-society/postdocs/abhilesh-dhawanjewar", "_blank");
        },
      },{
        id: 'social-x',
        title: 'X',
        section: 'Socials',
        handler: () => {
          window.open("https://twitter.com/abhilesh7", "_blank");
        },
      },{
      id: 'light-theme',
      title: 'Change theme to light',
      description: 'Change the theme of the site to Light',
      section: 'Theme',
      handler: () => {
        setThemeSetting("light");
      },
    },
    {
      id: 'dark-theme',
      title: 'Change theme to dark',
      description: 'Change the theme of the site to Dark',
      section: 'Theme',
      handler: () => {
        setThemeSetting("dark");
      },
    },
    {
      id: 'system-theme',
      title: 'Use system default theme',
      description: 'Change the theme of the site to System Default',
      section: 'Theme',
      handler: () => {
        setThemeSetting("system");
      },
    },];
