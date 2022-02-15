/// <reference types="react" />
import { ConnectedProps } from "react-redux";
export declare const LoadingIndicatorView: (props: any) => JSX.Element;
export declare const LoadingIndicatorDialog: (props: any) => JSX.Element;
declare const connector: import("react-redux").InferableComponentEnhancerWithProps<{
    setDataset: (value: any) => any;
}, {}>;
declare type Props = ConnectedProps<typeof connector>;
export declare const DatasetTabPanel: import("react-redux").ConnectedComponent<({ setDataset }: Props) => JSX.Element, Pick<{
    setDataset: (value: any) => any;
}, never>>;
export {};
//# sourceMappingURL=DatasetTabPanel.d.ts.map